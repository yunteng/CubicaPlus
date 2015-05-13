#include <util/SIMPLE_PARSER.h>
#include <util/MATRIX_UTIL.h>

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::PARTITIONED_HYBRID_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _fullDofs(tetMesh->partitionedFullsimDofs()),
  _reducedDofs(tetMesh->partitionedReducedsimDofs()),
  _rigger(rigger),
  _scd(NULL)
{
  //initialize the full space material force & matrix computation helper
  _fullMaterialCache = new FULL_MATERIAL_CACHE(tetMesh);

  //initialize the subspace material force & matrix computation helper
  _subMaterialCache = new SUB_MATERIAL_CACHE(tetMesh);

  /*
  functions in _fullMaterialCache are called
  before those in _subMaterialCache, let 
  _subMaterialCache be aware of to save redundant 
  function calls
  */ 
  _subMaterialCache->setCoupledWithFullspace(true);

  //initialize self collision detector
  _scd = new SELF_COLLISION_DETECTOR<BONE>(tetMesh, rigger);

  // read the self-collison penalty spring constant
  _scfMultiplier = SIMPLE_PARSER::getFloat("scf multiplier", 1.0);

  // read the domain interface spring constant
  _interfaceSpringConstant = SIMPLE_PARSER::getFloat("interface spring constant", 100.0);

  if(SIMPLE_PARSER::getBool("verbose", true)){
    cout << " SCF multiplier: " << _scfMultiplier << endl;
    cout << " Interface spring constant: " << _interfaceSpringConstant << endl;
  }

  // dynamic or quasistatic simulation?
  _dynamic = SIMPLE_PARSER::getBool("dynamic", false);

  /*
  resize a bunch of vectors to the number of skeletal partitions
  */
  _transformedUs.resize(_tetMesh->totalPartitions());
  _untransformedUs.resize(_tetMesh->totalPartitions());

  _UTKsf.resize(_tetMesh->totalPartitions());

  _previousFullDofs.resize(_tetMesh->totalPartitions(), 0);

  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    _transformedUs[x].resizeLike(_tetMesh->partitionBasis(x));
    _untransformedUs[x].resizeLike(_tetMesh->partitionBasis(x));
  }

  _reducedSCJacobians.resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());

  /*
  used by openmp to parallelize self collision 
  force jacobian computation 
  */
  _SCJacobianCopies = new MATRIX[_tetMesh->totalCores()];
  for(int x = 0; x < _tetMesh->totalCores(); x++){
    _SCJacobianCopies[x].resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  }

  // the subspace interface spring jacobian
  // is constant, so just comptue it once 
  // and store it
  computeReducedInterfaceSpringJacobians(_reducedInterfaceJacobians);


  if(_dynamic){
      _dt = SIMPLE_PARSER::getFloat("timestep", 1.0 / 60.0);
    if(SIMPLE_PARSER::getBool("verbose", true)){
      cout << " Dynamic simulation " << endl
           << "   timestep:      " << _dt << endl;
    }
    _tetMesh->recoverX();
    _positionOld = _tetMesh->x();
    _acceleration.resize(_tetMesh->dofs());
    _velocityOld.resize(_tetMesh->dofs());

    _acceleration.setZero();
    _velocityOld.setZero();

    _defaultOrderMass.resize(_tetMesh->dofs());

    /*
    the massVector in tetMesh store the mass per vertex, stretch it to all 3 dimensions
    */
    VECTOR& MVec = _tetMesh->massVector();
    for(int x = 0; x < MVec.size(); x++){
      _defaultOrderMass[x * 3] = MVec[x];
      _defaultOrderMass[x * 3 + 1] = MVec[x];
      _defaultOrderMass[x * 3 + 2] = MVec[x];
    }

    _reducedMasses.resize(_tetMesh->totalPartitions());
    _tetMesh->changeToPartitionOrder(_defaultOrderMass, _diagMass);

    // compute the reduced mass matrix for each partition
    for(int x = 0; x < _tetMesh->totalPartitions(); x++){
      int offset = _tetMesh->partitionDofStartIdx(x);
      int dofs = _tetMesh->partitionedVertices(x).size() * 3;
      _reducedMasses[x] = _tetMesh->partitionBasis(x).transpose() * _diagMass.segment(offset, dofs).asDiagonal() * _tetMesh->partitionBasis(x);
    }
    /*
    used by the dynamic oracles to keep track of 
    the vertices that should be kept in full space region
    */
    _keepDynamicRegion.resize(_tetMesh->unconstrainedNodes());
    _keepDynamicRegion.setZero();
  }
  /*
  used by the contact oracle, denote whether a vertex is in full space region
  */
  _isFullsim.conservativeResize(_tetMesh->unconstrainedNodes());
  _isFullsim.setZero();

  TIMING_BREAKDOWN::clear();
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::~PARTITIONED_HYBRID_INTEGRATOR()
{
  delete[] _SCJacobianCopies;

  if(_fullMaterialCache){
    delete _fullMaterialCache;
    _fullMaterialCache = NULL;
  }
  if(_subMaterialCache){
    delete _subMaterialCache;
    _subMaterialCache = NULL;
  }
  if(_scd){
    delete _scd;
    _scd = NULL;
  }
}
//////////////////////////////////////////
// compute the deformation gradients and 
// decompose them as see fit by the 
// material formulation
//////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeMaterialCache()
{
  _fullMaterialCache->cachePartialDecompositions();
  _subMaterialCache->cacheDecompositions();
}

//////////////////////////////////////////
// compute the full sim region around vertices in collision
//////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeFullsimRegions()
{
  // gather current in-collision vertices
  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();
  vector<pair<VEC3F*, SURFACE*> >& externalCollisionPairs = _tetMesh->collisionPairs();
  vector<EX_COLLISION_INFO>& externalCollisionPoints = _scd->externalCollisionPoints();

  Real maxDist = -1;

  TIMING_BREAKDOWN::tic();

  vector<int> contactVertices;

  for(unsigned int x = 0; x < externalCollisionPoints.size(); x++){
    VEC3I& triangleVertexIDs = externalCollisionPoints[x].triangleVertexIDs;
    for(int y = 0; y < 3; y++){
      contactVertices.push_back(triangleVertexIDs[y]);
    }
  }

  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){
    SELF_COLLISION_INFO& info = selfCollisionPoints[x];
    if(info.isCubaturePair){
      if(_isFullsim[info.vertexID] > 0 || _isFullsim[info.triangleVertexIDs[0]] || _isFullsim[info.triangleVertexIDs[1]] || _isFullsim[info.triangleVertexIDs[2]]){
        info.isCubaturePair = false;
        info.cubatureWeight = 1.0;
      }else{
        continue;
      }
    }

    contactVertices.push_back(info.vertexID);

    for(int y = 0; y < 3; y++){
      contactVertices.push_back(info.triangleVertexIDs[y]);
    }
  }

  for(unsigned int x = 0; x < externalCollisionPairs.size(); x++){
    int vertexID = _tetMesh->vertexID(externalCollisionPairs[x].first);

    contactVertices.push_back(vertexID);
  }

  vector<int> contactVerticesCopy = contactVertices;

  // gather in-collision vertices from the previous frame
  copy(_previousCollisionPoints.begin(), _previousCollisionPoints.end(), back_inserter(contactVertices));

  _previousCollisionPoints = contactVerticesCopy;

  TIMING_BREAKDOWN::toc("Gather Collision Vertices");


  TIMING_BREAKDOWN::tic();

  // compute full space region based on the in-collision vertices
  if(!_dynamic)
    _tetMesh->initPartitionedHybridSim(contactVertices, _isFullsim);
  else
    _tetMesh->initPartitionedHybridSim(contactVertices, _keepDynamicRegion, _isFullsim);

  TIMING_BREAKDOWN::toc("Init Partitioned Adaptive Mixed Sim");

  _subMaterialCache->registerFullsimTets();
}

//////////////////////////////////////////
// initialize a simulation step
//////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::initializeImplicitStep()
{
  static bool firstStep = true;
  if(firstStep){
    _positionOld = _tetMesh->x();
    firstStep = false;
  }

  TIMING_BREAKDOWN::tic();
  _tetMesh->recoverX();
  // store initial position in case we need to redo this frame
  _initPosition = _tetMesh->x();
  TIMING_BREAKDOWN::toc("Recover X");

  // cache the skinning rotations for the cubature tets
  TIMING_BREAKDOWN::tic();
  _subMaterialCache->cacheKeyTetTransforms(_rigger->skinningRotation());
  TIMING_BREAKDOWN::toc("Cache Key Tet Transforms");

  TIMING_BREAKDOWN::tic();
  _scd->vertexVsTetSCD();
  TIMING_BREAKDOWN::toc("Self Collision Detection");

  // compute full region using contact and dynamic oracle
  computeFullsimRegions();


  // check if any partition is entirely in full sim
  _hasEntireFullPartitions = false;
  for(unsigned int x = 0; x < _reducedDofs.size(); x++){
    if(_reducedDofs[x] == 0){
      _hasEntireFullPartitions = true;
      break;
    }
  }

  TIMING_BREAKDOWN::tic();


  // transform and reorder the basis for each partition based on the new ordering induced by condensation
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0)
      continue;

    int cols = _tetMesh->partitionBasis(x).cols();
    MATRIX& untransformedU = _tetMesh->partitionBasis(x);

    vector<int>& vertexIDs = _tetMesh->partitionedVertices(x);
    vector<int>& newOrderings = _tetMesh->adaptivePartitionedVertexOrdering(x);

    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(unsigned int y = 0; y < vertexIDs.size(); y++){
      int index = newOrderings[y];
      _transformedUs[x].block(index * 3, 0, 3, cols) = _rigger->skinningRotation()[vertexIDs[y]] * untransformedU.block(y * 3, 0, 3, cols);
      _untransformedUs[x].block(index * 3, 0, 3, cols) = untransformedU.block(y * 3, 0, 3, cols);
    }
  }
  TIMING_BREAKDOWN::toc("Compute transformed U");


  TIMING_BREAKDOWN::tic();
  VECTOR restDisp = _tetMesh->x() - _rigger->skinningDisp();
  vector<MATRIX3>& skinningRotation = _rigger->skinningRotation();


  // sync the subspace coordinates q with the newly updated (by skinning) full space displacement 
  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(unsigned int x = 0; x < skinningRotation.size(); x++){
    restDisp.segment<3>(x * 3) = skinningRotation[x].transpose() * restDisp.segment<3>(x * 3);
  }
  vector<VECTOR> perPartitionRestDisp;
  _tetMesh->changeToPartitionOrder(restDisp, perPartitionRestDisp);

  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0)
      _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) = _tetMesh->partitionBasis(x).transpose() * perPartitionRestDisp[x];
    else if(_reducedDofs[x] > 0)
      _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) = _untransformedUs[x].bottomRows(_reducedDofs[x]).transpose() * perPartitionRestDisp[x].tail(_reducedDofs[x]);
  }

  TIMING_BREAKDOWN::toc("Recover q");

  // change the orginal mesh displacement vector to partition order
  TIMING_BREAKDOWN::tic();
  _tetMesh->changeToPartitionOrder(_tetMesh->x(), _partitionedX);

  /*
  Some part of the current subspace region 
  might have been in full sim in the 
  previous frame so it may endure out-of-basis
  deformations, which cannot be properly removed 
  by the basis.
  We can stomp them away by recomputing the 
  displacement of the subspace region of each 
  partition using their subspace coordinates and
  skinning displacement.
  */
  _tetMesh->changeToPartitionOrder(_rigger->skinningDisp(), _workspace);

  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] != 0){
      _partitionedX.segment(_tetMesh->partitionDofStartIdx(x) + _fullDofs[x], _reducedDofs[x]) = _workspace.segment(_tetMesh->partitionDofStartIdx(x) + _fullDofs[x], _reducedDofs[x]) + _transformedUs[x].bottomRows(_reducedDofs[x]) * _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));
    }else{
      VECTOR dispDelta = _tetMesh->partitionBasis(x) * _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));
      vector<int>& vertexIDs = _tetMesh->partitionedVertices(x);

      #if USING_SUBSPACE_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for(unsigned int y = 0; y < vertexIDs.size(); y++){
        dispDelta.segment<3>(y * 3) = _rigger->skinningRotation()[vertexIDs[y]] * dispDelta.segment<3>(y * 3);
      }
      _partitionedX.segment(_tetMesh->partitionDofStartIdx(x), _reducedDofs[x]) = _workspace.segment(_tetMesh->partitionDofStartIdx(x), _reducedDofs[x]) + dispDelta;
    }
  }
  _tetMesh->restoreDefaultOrder(_partitionedX, _tetMesh->x());

  TIMING_BREAKDOWN::toc("Update Reduced Sim Positions");

  /*
  the interface spring jacobians (both full and 
  reduced) remains constant though a simulation
  frame. So we just need to compute them once.
  Add the mass matrix if we are doing dynamic simulation
  */
  TIMING_BREAKDOWN::tic();
  computeFixedSystemMatrices();
  TIMING_BREAKDOWN::toc("Compute Fixed System Matrices");

  if(_dynamic)
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  
}
////////////////////////////////////////////////
// the interface spring jacobians (both full and 
// reduced) remains constant though a simulation
// frame. So we just need to compute them once.
// Add the mass matrix if we are doing dynamic simulation
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeFixedSystemMatrices()
{
  // resize a bunch of matrices based on the 
  // current condensation partitioning
  BLOCK_SPARSE_MATRIX& partitionedKff = _tetMesh->partitionedFullStiffnessDiag();
  partitionedKff.resizeAndWipe(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());
  partitionedKff.setBlockDimensions(_fullDofs, _fullDofs);

  BLOCK_SPARSE_MATRIX& partitionedKsf = _tetMesh->partitionedFullStiffnessOffDiag();
  partitionedKsf.resizeAndWipe(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());
  partitionedKsf.setBlockDimensions(_reducedDofs, _fullDofs);

  BLOCK_SPARSE_MATRIX& partitionedBoundaryKss = _tetMesh->partitionedReducedStiffnessDiag();
  partitionedBoundaryKss.resizeAndWipe(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());
  partitionedBoundaryKss.setBlockDimensions(_reducedDofs, _reducedDofs);


  BLOCK_COO_MATRIX fixedFullHessian(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());
  fixedFullHessian.setBlockDimensions(_fullDofs, _fullDofs);

  _fixedReducedHessian.resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  _fixedReducedHessian.setZero();

  _SCJacobians.resizeAndWipe(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());
  _SCJacobians.setBlockDimensions(_fullDofs, _fullDofs);

  // compute the interface spring jacobians
  computeInterfaceSpringJacobians(fixedFullHessian, _fixedReducedHessian);

  if(_dynamic){

    _tetMesh->changeToPartitionOrder(_defaultOrderMass, _diagMass);

    // add the mass matrix
    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int x = 0; x < _tetMesh->totalPartitions(); x++){
      
      int reducedOffset = _tetMesh->partitionRankStartIdx(x);
      int rank = _tetMesh->partitionRank(x);

      if(_fullDofs[x] == 0){
        _fixedReducedHessian.block(reducedOffset, reducedOffset, rank, rank) += (1.0 / _dt / _dt) * _reducedMasses[x];
        continue;
      }

      int fullOffset = _tetMesh->partitionDofStartIdx(x);
      COO_MATRIX& mat = fixedFullHessian(x, x);
      for(int y = 0; y < _fullDofs[x]; y++){
        mat.add(_diagMass[fullOffset + y] / _dt / _dt, y, y);
      }

      MATRIX reducedMass = _transformedUs[x].bottomRows(_reducedDofs[x]).transpose() * _diagMass.segment(fullOffset + _fullDofs[x], _reducedDofs[x]).asDiagonal() * _transformedUs[x].bottomRows(_reducedDofs[x]);

      _fixedReducedHessian.block(reducedOffset, reducedOffset, rank, rank) += (1.0 / _dt / _dt) * reducedMass;
    }
  }

  // convert to Eigen format for fast matrix-vector multiplication
  COO_MATRIX tmp;
  fixedFullHessian.toCOOMatrix(tmp);
  tmp.toSpMat(_fixedFullHessian);

  _stiffnessDiagonal = _fixedFullHessian.diagonal();
}
////////////////////////////////////////////////
// finalize the simualtion step, update the mesh 
// if we are doing dynamic simulation, register 
// the vertices that need to be kept in full sim
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::finalizeImplicitStep()
{
  _tetMesh->updateFullMesh();

  _previousFullDofs  = _fullDofs;


  // dynamic oracle, check if the velocities or 
  // accelerations of the current full sim regions 
  // are small enough so that we can deactivate them
  if(_dynamic){
    _velocityOld = (1.0 / _dt) * (_tetMesh->x() - _positionOld);
    _positionOld = _tetMesh->x();
    
    _keepDynamicRegion.setZero();

    if(activeFullDofs() > 0){

      VECTOR partitionedAcceleration;

      _tetMesh->changeToPartitionOrder(_velocityOld, _partitionedVelocity);

      _tetMesh->changeToPartitionOrder(_acceleration, partitionedAcceleration);

      for(int x = 0; x < _tetMesh->totalPartitions(); x++){
        if(_fullDofs[x] == 0)
          continue;

        VECTOR velocity = _partitionedVelocity.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x] + _reducedDofs[x]);

        Real fullsimAvgVelocity = 0;
        for(int y = 0; y < _fullDofs[x] / 3; y++){
          fullsimAvgVelocity += velocity.segment<3>(y * 3).norm();
        }
        fullsimAvgVelocity /= _fullDofs[x] / 3;

        Real reducedsimAvgVelocity = 0;
        for(int y = _fullDofs[x] / 3; y < velocity.size() / 3; y++){
          reducedsimAvgVelocity += velocity.segment<3>(y * 3).norm();
        }
        reducedsimAvgVelocity /= _reducedDofs[x] / 3;

        Real diff = abs(fullsimAvgVelocity - reducedsimAvgVelocity) / fullsimAvgVelocity;

        VECTOR accleration = partitionedAcceleration.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x] + _reducedDofs[x]);

        Real fullsimAvgAcceleration = 0;
        for(int y = 0; y < _fullDofs[x] / 3; y++){
          fullsimAvgAcceleration += accleration.segment<3>(y * 3).norm();
        }
        fullsimAvgAcceleration /= _fullDofs[x] / 3;

        Real reducedsimAvgAcceleration = 0;
        for(int y = _fullDofs[x] / 3; y < velocity.size() / 3; y++){
          reducedsimAvgAcceleration += accleration.segment<3>(y * 3).norm();
        }
        reducedsimAvgAcceleration /= _reducedDofs[x] / 3;

        Real accDiff = abs(fullsimAvgAcceleration - reducedsimAvgAcceleration) / fullsimAvgAcceleration;

        if(fullsimAvgAcceleration > 0.2 || fullsimAvgAcceleration > 0.2){
          vector<int>& vertexIDs = _tetMesh->partitionedVertices(x);
          vector<int>& newOrderings = _tetMesh->adaptivePartitionedVertexOrdering(x);
          int fullsimVertices = _fullDofs[x] / 3;
          for(unsigned int y = 0; y < vertexIDs.size(); y++){
            if(newOrderings[y] < fullsimVertices)
              _keepDynamicRegion[vertexIDs[y]] = 1.0;
          }
        } 
      } 
    }
  }
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::setPosition(VECTOR& newPosition)
{
  // update the positions for the full sim regions
  if(activeFullDofs() > 0){
    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int x = 0; x < _tetMesh->totalPartitions(); x++){
      if(_fullDofs[x] == 0)
        continue;
      _partitionedX.segment(_tetMesh->partitionDofStartIdx(x), _fullDofs[x]) -= newPosition.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]);
    }
    _qi = newPosition.tail(_reducedRegionGradient.size());
  }

  /*
  compute the displacement delta of the subspace region using _qi
  */
  vector<VECTOR> dispDeltas(_tetMesh->totalPartitions());

  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0){
      dispDeltas[x] = _tetMesh->partitionBasis(x) * _qi.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));

      vector<int>& vertexIDs = _tetMesh->partitionedVertices(x);
      
      for(unsigned int y = 0; y < vertexIDs.size(); y++){
        dispDeltas[x].segment<3>(y * 3) = _rigger->skinningRotation()[vertexIDs[y]] * dispDeltas[x].segment<3>(y * 3);
      }
    }
  }
  // update the subspace region
  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0){
      _partitionedX.segment(_tetMesh->partitionDofStartIdx(x), _reducedDofs[x]) -= dispDeltas[x];
    }else{
      _partitionedX.segment(_tetMesh->partitionDofStartIdx(x) + _fullDofs[x], _reducedDofs[x]) -= _transformedUs[x].bottomRows(_reducedDofs[x]) * _qi.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));
    }
  }

  _tetMesh->q() -= _qi;

  // map back to ogrinal mesh displacement vector
  _tetMesh->restoreDefaultOrder(_partitionedX, _tetMesh->x());

  if(_dynamic){
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }
}
////////////////////////////////////////////////
// compute the energy of the system
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
Real PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeSystemEnergy()
{
  /*
  approximate the elastic energy using the 
  material force cubatures. This may not be very 
  accurate if we have full space regions, but 
  currently we are not using the energy directly 
  for the Newton process
  */
  _energy = _subMaterialCache->computeElasticEnergy();

  // add kenematic energy
  if(_dynamic)
    _energy += 0.5 * _dt * _dt * _acceleration.dot(_defaultOrderMass.asDiagonal() * _acceleration);
  
  return _energy;
}
////////////////////////////////////////////////
// compute the system matrices, partitioned Kff, 
// Ksf and Kss
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
BLOCK_SPARSE_MATRIX& PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeSystemMatrix()
{
  BLOCK_SPARSE_MATRIX& partitionedKff = _tetMesh->partitionedFullStiffnessDiag();
  BLOCK_SPARSE_MATRIX& partitionedKsf = _tetMesh->partitionedFullStiffnessOffDiag();
  BLOCK_SPARSE_MATRIX& partitionedBoundaryKss = _tetMesh->partitionedReducedStiffnessDiag();

  partitionedKff.clear();
  partitionedKsf.clear();
  partitionedBoundaryKss.clear();

  TIMING_BREAKDOWN::tic();
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0)
      continue;
    COO_MATRIX Kff;
    COO_MATRIX Ksf;
    COO_MATRIX BoundaryKss;

    _fullMaterialCache->computePartialStiffnessMatrix(x, Kff, Ksf, BoundaryKss);
    partitionedKff.equals(Kff, x, x);
    partitionedKsf.equals(Ksf, x, x);
    partitionedBoundaryKss.equals(BoundaryKss, x, x);

    // U_s^T * K_sf
    _UTKsf[x] = (_transformedUs[x].bottomRows(_reducedDofs[x]).transpose() * partitionedKsf(x, x)).sparseView();

    _stiffnessDiagonal.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]) += partitionedKff(x, x).diagonal();
  }
  TIMING_BREAKDOWN::toc("Compute Full Stiffness Matrix");


  TIMING_BREAKDOWN::tic();

  _reducedKss = _subMaterialCache->computeReducedStiffnessMatrix();


  TIMING_BREAKDOWN::toc("Compute Reduced Stiffness Matrix");

  TIMING_BREAKDOWN::tic();

  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0)
      continue;

    int offSet = _tetMesh->partitionRankStartIdx(x);
    int rank = _tetMesh->partitionRank(x);

    _reducedKss.block(offSet, offSet, rank, rank) += _transformedUs[x].bottomRows(_reducedDofs[x]).transpose() * partitionedBoundaryKss(x, x) * _transformedUs[x].bottomRows(_reducedDofs[x]);
  }
  TIMING_BREAKDOWN::toc("Correct Reduced Stiffness Matrix");

  TIMING_BREAKDOWN::tic();
  _reducedKss += _fixedReducedHessian;
  TIMING_BREAKDOWN::toc("Add Fixed Internal Hessian");

  computeCollisionMatrices();

  TIMING_BREAKDOWN::tic();
  // if any partition is entirely in full sim
  // _reducedKss would have zero diagonal blocks,
  // and the matrix inversion would fail
  // remove those blocks by shrinking the matrix
  // if necessary
  if(_hasEntireFullPartitions){
    removeZeroBlocks(_reducedKss, _prunedReducedKss);
    _reducedKssInv.compute(_prunedReducedKss);
  }else{
    _reducedKssInv.compute(_reducedKss);
  }
  TIMING_BREAKDOWN::toc("Invert Reduced Internal Hessian");

  return partitionedKff;
}

////////////////////////////////////////////////
// compute collision matrices
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeCollisionMatrices()
{
  TIMING_BREAKDOWN::tic();
  // self collison jacobian
  computeSelfCollisionSpringForceJacobian(_SCJacobians, _reducedKss);
  // external collision jacobian
  computeExternalCollisionForceJacobian(_SCJacobians);

  // convert to eigen format for fat matrix-vector
  // multiplication
  COO_MATRIX tmp;
  _SCJacobians.toCOOMatrix(tmp);
  tmp.toSpMat(_SCJacobiansSpMat);

  // update the diagonal of the entire system matrix
  _hessianDiagonal = _stiffnessDiagonal + _SCJacobiansSpMat.diagonal();

  // add the interface spring jocobians
  // (and mass matrix if dynamic)
  _SCJacobiansSpMat += _fixedFullHessian;
  _SCJacobiansSpMat.makeCompressed();

  TIMING_BREAKDOWN::toc("Compute Collision Jacobians");
}
////////////////////////////////////////////////
// This is called when we are doing subspace-only 
// simulation, compute qi using prefactored 
// system matrix
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeQi(const VECTOR& u_s)
{
  _reducedKssInv.solve(u_s, _qi);
}

////////////////////////////////////////////////
// Compute the RHS of the linear system
////////////////////////////////////////////////
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
VECTOR& PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeSystemForce()
{
  // material force of the full space region
  vector<VECTOR>& fullForces = _fullMaterialCache->computePartialInternalForce();

  _fullRegionGradient.conservativeResize(_tetMesh->partitionFullsimDofStartIdx(_tetMesh->totalPartitions()));

  // reduced material force of the subspace region
  _reducedRegionGradient = _subMaterialCache->computeInternalForce();
  _reducedRegionGradient *= -1.0;

  // inertia force
  if(_dynamic){ 
    _tetMesh->changeToPartitionOrder(_defaultOrderMass.asDiagonal() * _acceleration, _accelerationForces);
  }  

  // subspace boundary force
  TIMING_BREAKDOWN::tic();
  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] > 0){
      _reducedRegionGradient.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) += _transformedUs[x].bottomRows(_reducedDofs[x]).transpose() * fullForces[x].tail(_reducedDofs[x]);

      _fullRegionGradient.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]) = fullForces[x].head(_fullDofs[x]); 
    }
    if(_dynamic){
      if(_fullDofs[x] == 0){
        vector<int>& vertexIDs = _tetMesh->partitionedVertices(x);
        VECTOR& accForce = _accelerationForces[x];
        for(unsigned int y = 0; y < vertexIDs.size(); y++){
          accForce.segment<3>(y * 3) = _rigger->skinningRotation()[vertexIDs[y]].transpose() * accForce.segment<3>(y * 3);
        }
        _reducedRegionGradient.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) += _tetMesh->partitionBasis(x).transpose() * accForce;

      }else{
        _reducedRegionGradient.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) += _transformedUs[x].bottomRows(_reducedDofs[x]).transpose() * _accelerationForces[x].tail(_reducedDofs[x]);

        _fullRegionGradient.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]) += _accelerationForces[x].head(_fullDofs[x]);
      }
      
    }
  }
  TIMING_BREAKDOWN::toc("projection force");

  // self collision force
  TIMING_BREAKDOWN::tic();
  _energy += computeSelfCollisionSpringForces(_fullRegionGradient, _reducedRegionGradient);

  // external collision force
  _energy += computeExternalCollisionForces(_fullRegionGradient);

  TIMING_BREAKDOWN::toc("Compute Collision Forces");

  TIMING_BREAKDOWN::tic();
  // interface spring force
  _energy += computeInterfaceSpringForces(_fullRegionGradient, _reducedRegionGradient);

  TIMING_BREAKDOWN::toc("Compute Coupling Forces");

  // pack everything into one big vector
  _gradient.conservativeResize(_fullRegionGradient.size() + _reducedRegionGradient.size());
  _gradient.head(_fullRegionGradient.size()) = _fullRegionGradient;
  _gradient.tail(_reducedRegionGradient.size()) = _reducedRegionGradient;
  
  return _gradient;
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
Real PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForces(VECTOR& fullForceVector, VECTOR& reducedForceVector)
{
  Real energy = 0;

  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();

  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){
    SELF_COLLISION_INFO& info = selfCollisionPoints[x];

    int leftIndex = info.vertexID;

    VEC3F& leftVertex = _tetMesh->vertices()[leftIndex];
    VEC3I& triangleIndices = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    Real multiplier = _scfMultiplier * info.avgArea * info.cubatureWeight;

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndices[x]));
    }

    MATRIX3& M = info.M;

    VEC3F penaltyForce = info.M * (surfacePosition - leftVertex) * multiplier;

    /*
    TODO::add dynamic damping
    */


    energy += 0.5 * (surfacePosition - leftVertex).dot(penaltyForce);

    pair<int, int>& partitionID = _tetMesh->partitionedVertexID(leftIndex);

    if(!info.isCubaturePair){
      int newID = _tetMesh->adaptivePartitionedVertexID(partitionID);
        fullForceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(partitionID.first) + newID * 3) -= penaltyForce;

      for(int i = 0; i < 3; i++){
        pair<int, int>& partitionID = _tetMesh->partitionedVertexID(triangleIndices[i]);

        int newID = _tetMesh->adaptivePartitionedVertexID(partitionID);

        fullForceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(partitionID.first) + newID * 3) += lambda[i] * penaltyForce;
      }
    }else{
      reducedForceVector.segment(_tetMesh->partitionRankStartIdx(partitionID.first), _tetMesh->partitionRank(partitionID.first))
       -= _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * (_rigger->skinningRotation()[leftIndex].transpose() * penaltyForce);

     for(int i = 0; i < 3; i++){
        pair<int, int>& partitionID = _tetMesh->partitionedVertexID(triangleIndices[i]);

        reducedForceVector.segment(_tetMesh->partitionRankStartIdx(partitionID.first), _tetMesh->partitionRank(partitionID.first))
           += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * (_rigger->skinningRotation()[triangleIndices[i]].transpose() * (lambda[i] * penaltyForce));
      }
    }
  }

  return energy;
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
Real PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeExternalCollisionForces(VECTOR& forceVector)
{
  Real energy = 0;
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();
  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    pair<int, int>& partitionID = _tetMesh->partitionedVertexID(vertexID);
    
    int newID = _tetMesh->adaptivePartitionedVertexID(partitionID);

    VEC3F force = surface->force(*vertex);
    forceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(partitionID.first) + newID * 3) += force;
    energy += 0.5 * force.squaredNorm() / surface->collisionStiffness();
  }

  vector<EX_COLLISION_INFO>& externalCollisionPoints = _scd->externalCollisionPoints();

  for(unsigned int x = 0; x < externalCollisionPoints.size(); x++){

    EX_COLLISION_INFO& info = externalCollisionPoints[x];

    VEC3F& leftVertex = info.penetratingPosition;

    VEC3I& triangleIndices = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    Real multiplier = _scfMultiplier * info.avgArea;

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndices[x]));
    } 

    MATRIX3& M = info.M;

    VEC3F penaltyForce = info.M * (surfacePosition - leftVertex) * multiplier;

    energy += 0.5 * (surfacePosition - leftVertex).dot(penaltyForce);

    for(int i = 0; i < 3; i++){
      pair<int, int>& partitionID = _tetMesh->partitionedVertexID(triangleIndices[i]);

      int newID = _tetMesh->adaptivePartitionedVertexID(partitionID);

      forceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(partitionID.first) + newID * 3) += lambda[i] * penaltyForce;
    }
  }
  return energy;
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForceJacobian(BLOCK_COO_MATRIX& fullSystemMatrix, MATRIX& reducedSystemMatrix)
{
  fullSystemMatrix.clear();

  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();

  int cubatureCnt = 0;
  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){

    SELF_COLLISION_INFO& info = selfCollisionPoints[x];

    if(info.isCubaturePair){
      cubatureCnt++;
      continue;
    }  

    int leftIndex = info.vertexID;

    VEC3F& leftVertex = _tetMesh->vertices()[leftIndex];
    VEC3I& triangleIndices = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    Real multiplier = _scfMultiplier * info.avgArea;

    // if(_dynamic)
      // multiplier += _scDampingConst * _alpha[3];

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndices[x]));
    } 

    VEC3F normal = surfacePosition - leftVertex;

    normal.normalize();

    MATRIX3 M = MATRIX_UTIL::outer_product(normal);
    M *= 0.5;
    M += 0.5 * MATRIX3::Identity();
    info.M = M;

    MATRIX3 leftJocobian = multiplier * M;

    pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(leftIndex);

    COO_MATRIX& mat = fullSystemMatrix(leftPartitionID.first, leftPartitionID.first);

    int leftNewID = _tetMesh->adaptivePartitionedVertexID(leftPartitionID);

    mat.add3x3(leftJocobian, leftNewID * 3, leftNewID * 3);

    for(int x = 0; x < 3; x++){

      pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);

      int rightNewID = _tetMesh->adaptivePartitionedVertexID(rightPartitionID);

      MATRIX jacobian = -lambda[x] * multiplier * M;
      
      COO_MATRIX& Bxy = fullSystemMatrix(leftPartitionID.first, rightPartitionID.first);

      Bxy.add3x3(jacobian, leftNewID * 3, rightNewID * 3);

      COO_MATRIX& Byx = fullSystemMatrix(rightPartitionID.first, leftPartitionID.first);

      Byx.add3x3(jacobian.transpose(), rightNewID * 3, leftNewID * 3);
    }

    for(int x = 0; x < 3; x++){
      pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);

      int leftNewID = _tetMesh->adaptivePartitionedVertexID(leftPartitionID);

      for(int y = 0; y < 3; y++){
        pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[y]);
 
        MATRIX3 rightJacobian = lambda[x] * lambda[y] * multiplier * M;

        COO_MATRIX& mat = fullSystemMatrix(leftPartitionID.first, rightPartitionID.first);

        int rightNewID = _tetMesh->adaptivePartitionedVertexID(rightPartitionID);

        mat.add3x3(rightJacobian, leftNewID * 3, rightNewID * 3);
      }
    }
  }

  if(cubatureCnt == 0)
    return;

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_SUBSPACE_OPENMP
    const int id  = omp_get_thread_num();
    #else
    const int id = 0;
    #endif
    _SCJacobianCopies[id].setZero();

    #if USING_SUBSPACE_OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){

      SELF_COLLISION_INFO& info = selfCollisionPoints[x];

      if(!info.isCubaturePair)
        continue;

      int leftIndex = info.vertexID;

      VEC3F& leftVertex = _tetMesh->vertices()[leftIndex];
      VEC3I& triangleIndices = info.triangleVertexIDs;
      VEC3F& lambda = info.baryCenter;

      Real multiplier = _scfMultiplier * info.avgArea * info.cubatureWeight;

      // if(_dynamic)
        // multiplier += _scDampingConst * _alpha[3];

      VEC3F surfacePosition;
      surfacePosition.setZero();
      for(int x = 0; x < 3; x++){
        surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndices[x]));
      } 

      VEC3F normal = surfacePosition - leftVertex;

      normal.normalize();

      MATRIX3 M = MATRIX_UTIL::outer_product(normal);
      M *= 0.5;
      M += 0.5 * MATRIX3::Identity();
      info.M = M;

      MATRIX3 leftJocobian = multiplier * M;

      pair<int, int>& partitionID = _tetMesh->partitionedVertexID(leftIndex);

      leftJocobian = _rigger->skinningRotation()[leftIndex].transpose() * leftJocobian * _rigger->skinningRotation()[leftIndex];

      int rowStart = _tetMesh->partitionRankStartIdx(partitionID.first);
      int colStart = rowStart;
      int subRank = _tetMesh->partitionRank(partitionID.first);

      _SCJacobianCopies[id].block(rowStart, colStart, subRank, subRank) += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * leftJocobian * _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second);

      for(int x = 0; x < 3; x++){
        pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);
        int rowStart = _tetMesh->partitionRankStartIdx(leftPartitionID.first);
        int rowRank  = _tetMesh->partitionRank(leftPartitionID.first);

        for(int y = 0; y < 3; y++){
          pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[y]);
          int colStart = _tetMesh->partitionRankStartIdx(rightPartitionID.first);
          int colRank  = _tetMesh->partitionRank(rightPartitionID.first);

          MATRIX3 rightJacobian = lambda[x] * lambda[y] * multiplier * M;

          rightJacobian = _rigger->skinningRotation()[triangleIndices[x]].transpose() * rightJacobian * _rigger->skinningRotation()[triangleIndices[y]];

          _SCJacobianCopies[id].block(rowStart, colStart, rowRank, colRank) += 
            _tetMesh->partitionVertexBasis(leftPartitionID.first, leftPartitionID.second).transpose() * rightJacobian * _tetMesh->partitionVertexBasis(rightPartitionID.first, rightPartitionID.second);
        }
      }

      pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(leftIndex);
      rowStart = _tetMesh->partitionRankStartIdx(leftPartitionID.first);
      int rowRank = _tetMesh->partitionRank(leftPartitionID.first);

      for(int x = 0; x < 3; x++){

        pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);

        int colStart = _tetMesh->partitionRankStartIdx(rightPartitionID.first);
        int colRank  = _tetMesh->partitionRank(rightPartitionID.first);

        MATRIX jacobian = -lambda[x] * multiplier * M;
        
        jacobian = _rigger->skinningRotation()[leftIndex].transpose() * jacobian * _rigger->skinningRotation()[triangleIndices[x]];

        _SCJacobianCopies[id].block(rowStart, colStart, rowRank, colRank) += _tetMesh->partitionVertexBasis(leftPartitionID.first, leftPartitionID.second).transpose() * jacobian * _tetMesh->partitionVertexBasis(rightPartitionID.first, rightPartitionID.second);

        _SCJacobianCopies[id].block(colStart, rowStart, colRank, rowRank) += _tetMesh->partitionVertexBasis(rightPartitionID.first, rightPartitionID.second).transpose() * jacobian.transpose() * _tetMesh->partitionVertexBasis(leftPartitionID.first, leftPartitionID.second);
      }
    }
  }
  for(int x = 0; x < _tetMesh->totalCores(); x++)
    reducedSystemMatrix += _SCJacobianCopies[x];
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeExternalCollisionForceJacobian(BLOCK_COO_MATRIX& systemMatrix)
{
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();

  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    pair<int, int>& partitionID = _tetMesh->partitionedVertexID(vertexID);
    
    int newID = _tetMesh->adaptivePartitionedVertexID(partitionID);

    MATRIX3 springJacobian = surface->springJacobian(*vertex);
    
    COO_MATRIX& mat = systemMatrix(partitionID.first, partitionID.first);

    mat.add3x3(springJacobian, newID * 3, newID * 3);
  }

  vector<EX_COLLISION_INFO>& externalCollisionPoints = _scd->externalCollisionPoints();

  for(unsigned int x = 0; x < externalCollisionPoints.size(); x++){

    EX_COLLISION_INFO& info = externalCollisionPoints[x];

    VEC3F& leftVertex = info.penetratingPosition;

    VEC3I& triangleIndices = info.triangleVertexIDs;
    VEC3F& lambda = info.baryCenter;

    Real multiplier = _scfMultiplier * info.avgArea;

    VEC3F surfacePosition;
    surfacePosition.setZero();
    for(int x = 0; x < 3; x++){
      surfacePosition += lambda[x] * *(_tetMesh->vertex(triangleIndices[x]));
    } 

    VEC3F normal = surfacePosition - leftVertex;

    normal.normalize();

    MATRIX3 M = MATRIX_UTIL::outer_product(normal);
    M *= 0.5;
    M += 0.5 * MATRIX3::Identity();
    info.M = M;

    for(int x = 0; x < 3; x++){
      pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);

      int leftNewID = _tetMesh->adaptivePartitionedVertexID(leftPartitionID);

      for(int y = 0; y < 3; y++){
        pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[y]);
 
        MATRIX3 rightJacobian = lambda[x] * lambda[y] * multiplier * M;

        COO_MATRIX& mat = systemMatrix(leftPartitionID.first, rightPartitionID.first);

        int rightNewID = _tetMesh->adaptivePartitionedVertexID(rightPartitionID);

        mat.add3x3(rightJacobian, leftNewID * 3, rightNewID * 3);
      }
    }
  }
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::getMatVecMult(const VECTOR& input, VECTOR& output)
{
  VECTOR q = input.tail(_reducedRegionGradient.size());
  VECTOR u = input.head(_fullRegionGradient.size());

  output.conservativeResize(input.size());
  output.head(_fullRegionGradient.size()) = _SCJacobiansSpMat * u;
  output.tail(_reducedRegionGradient.size()) = _reducedKss * q;

  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_fullDofs[x] == 0)
      continue;

    output.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]) += _tetMesh->partitionedFullStiffnessDiag()(x, x) * u.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]);

    if(_reducedDofs[x] > 0){
      output.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]) += _UTKsf[x].transpose() * q.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));

      output.segment(_fullRegionGradient.size() + _tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) += _UTKsf[x] * u.segment(_tetMesh->partitionFullsimDofStartIdx(x), _fullDofs[x]);
    }
  }
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeInterfaceSpringJacobians(BLOCK_COO_MATRIX& fullSystemMatrix, MATRIX& reducedSystemMatrix)
{
  reducedSystemMatrix = _reducedInterfaceJacobians;

  map<pair<int, int>, vector<pair<int, int> > > interfaceSprings = _tetMesh->interfaceSprings();

  MATRIX3 springJocobian = _interfaceSpringConstant * MATRIX3::Identity();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = interfaceSprings.begin(); iter != interfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;

    if(_fullDofs[leftPartition] == 0 || _fullDofs[rightPartition] == 0)
      continue;

    vector<int>& leftAdaptiveOrdering = _tetMesh->adaptivePartitionedVertexOrdering(leftPartition);
    vector<int>& rightAdaptiveOrdering = _tetMesh->adaptivePartitionedVertexOrdering(rightPartition);

    COO_MATRIX& Bxx = fullSystemMatrix(leftPartition, leftPartition);

    COO_MATRIX& Bxy = fullSystemMatrix(leftPartition, rightPartition);

    COO_MATRIX& Byx = fullSystemMatrix(rightPartition, leftPartition);
    
    COO_MATRIX& Byy = fullSystemMatrix(rightPartition, rightPartition);

    vector<pair<int, int> >& springs = iter->second;
    for(unsigned int x = 0; x < springs.size(); x++){
      int leftVID = springs[x].first;
      int rightVID = springs[x].second;

      if(_tetMesh->isPartitionedFullsimVertex(leftPartition, leftVID)){

        int rowStart = _tetMesh->partitionRankStartIdx(leftPartition);
        int colStart = _tetMesh->partitionRankStartIdx(rightPartition);
        int rowRank = _tetMesh->partitionRank(leftPartition);
        int colRank = _tetMesh->partitionRank(rightPartition);

        MATRIX leftVertexU = _tetMesh->partitionVertexBasis(leftPartition, leftVID);
        MATRIX rightVertexU = _tetMesh->partitionVertexBasis(rightPartition, rightVID);

        reducedSystemMatrix.block(rowStart, rowStart, rowRank, rowRank) -= leftVertexU.transpose() * springJocobian * leftVertexU;

        reducedSystemMatrix.block(colStart, colStart, colRank, colRank) -= rightVertexU.transpose() * springJocobian * rightVertexU;

        reducedSystemMatrix.block(rowStart, colStart, rowRank, colRank) += leftVertexU.transpose() * springJocobian * rightVertexU;

        reducedSystemMatrix.block(colStart, rowStart, colRank, rowRank) += rightVertexU.transpose() * springJocobian * leftVertexU;


        leftVID = leftAdaptiveOrdering[leftVID];
        rightVID = rightAdaptiveOrdering[rightVID];

        Bxx.add3x3(springJocobian, leftVID * 3, leftVID * 3);
        Byx.add3x3(-springJocobian, rightVID * 3, leftVID * 3);

        Byy.add3x3(springJocobian, rightVID * 3, rightVID * 3);
        Bxy.add3x3(-springJocobian, leftVID * 3, rightVID * 3);
      }
    }
  }
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
Real PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeInterfaceSpringForces(VECTOR& fullForceVector, VECTOR& reducedForceVector)
{
  Real energy = 0;
  
  map<pair<int, int>, vector<pair<int, int> > >& interfaceSprings = _tetMesh->interfaceSprings();

  map<pair<int, int>, vector<MATRIX> >& interfaceUs = _tetMesh->interfaceUs();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = interfaceSprings.begin(); iter != interfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;

    if(_fullDofs[leftPartition] == 0 && _fullDofs[rightPartition] == 0){
      vector<MATRIX>& springUs = interfaceUs[iter->first];
      VECTOR leftq = _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(leftPartition), _tetMesh->partitionRank(leftPartition));
      VECTOR rightq = _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(rightPartition), _tetMesh->partitionRank(rightPartition));

      VECTOR leftDisp = springUs[0] * leftq;
      VECTOR rightDisp = springUs[1] * rightq;
      VECTOR diff = leftDisp - rightDisp;

      energy += diff.squaredNorm() * 0.5 * _interfaceSpringConstant;

      leftq *= _interfaceSpringConstant;
      rightq *= _interfaceSpringConstant;

      reducedForceVector.segment(_tetMesh->partitionRankStartIdx(leftPartition), _tetMesh->partitionRank(leftPartition)) += springUs[2] * leftq - springUs[3] * rightq;

      reducedForceVector.segment(_tetMesh->partitionRankStartIdx(rightPartition), _tetMesh->partitionRank(rightPartition)) -= springUs[3].transpose() * leftq - springUs[4] * rightq;
      continue;
    }

    vector<int>& leftAdaptiveOrdering = _tetMesh->adaptivePartitionedVertexOrdering(leftPartition);
    vector<int>& rightAdaptiveOrdering = _tetMesh->adaptivePartitionedVertexOrdering(rightPartition);

    vector<pair<int, int> >& springs = iter->second;
    for(unsigned int x = 0; x < springs.size(); x++){
      int originalLeftVID = springs[x].first;
      int originalRightVID = springs[x].second;

      int leftVID = leftAdaptiveOrdering[originalLeftVID];
      int rightVID = rightAdaptiveOrdering[originalRightVID];

      VEC3F leftDisp = _partitionedX.segment<3>(_tetMesh->partitionDofStartIdx(leftPartition) + leftVID * 3);
      VEC3F rightDisp = _partitionedX.segment<3>(_tetMesh->partitionDofStartIdx(rightPartition) + rightVID * 3);

      VEC3F diff = _interfaceSpringConstant * (leftDisp - rightDisp);

      Real e = 0.5 * diff.dot(leftDisp - rightDisp);
      energy += e;

      if(e < 1e-6)
        continue;

      if(leftVID * 3 >= _fullDofs[leftPartition]){
        int vertexID = _tetMesh->partitionedVertices(leftPartition)[originalLeftVID];
        VEC3F force = _rigger->skinningRotation()[vertexID].transpose() * diff;

        reducedForceVector.segment(_tetMesh->partitionRankStartIdx(leftPartition), _tetMesh->partitionRank(leftPartition)) += _tetMesh->partitionVertexBasis(leftPartition, originalLeftVID).transpose() * force;
        reducedForceVector.segment(_tetMesh->partitionRankStartIdx(rightPartition), _tetMesh->partitionRank(rightPartition)) -= _tetMesh->partitionVertexBasis(rightPartition, originalRightVID).transpose() * force;
        continue;
      }
        

      fullForceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(leftPartition) + leftVID * 3) += diff;

      fullForceVector.segment<3>(_tetMesh->partitionFullsimDofStartIdx(rightPartition) + rightVID * 3) -= diff;      
    }
  }
  cout << "Interface spring energy " << energy << endl;

  return energy;
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::computeReducedInterfaceSpringJacobians(MATRIX& reducedSystemMatrix)
{
  reducedSystemMatrix.conservativeResize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  reducedSystemMatrix.setZero();
  map<pair<int, int>, vector<pair<int, int> > > interfaceSprings = _tetMesh->interfaceSprings();

  MATRIX3 springJocobian = _interfaceSpringConstant * MATRIX3::Identity();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = interfaceSprings.begin(); iter != interfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;

    vector<pair<int, int> >& springs = iter->second;
    for(unsigned int x = 0; x < springs.size(); x++){
      int leftVID = springs[x].first;
      int rightVID = springs[x].second;
      
      int rowStart = _tetMesh->partitionRankStartIdx(leftPartition);
      int colStart = _tetMesh->partitionRankStartIdx(rightPartition);
      int rowRank = _tetMesh->partitionRank(leftPartition);
      int colRank = _tetMesh->partitionRank(rightPartition);

      MATRIX leftVertexU = _tetMesh->partitionVertexBasis(leftPartition, leftVID);
      MATRIX rightVertexU = _tetMesh->partitionVertexBasis(rightPartition, rightVID);

      reducedSystemMatrix.block(rowStart, rowStart, rowRank, rowRank) += leftVertexU.transpose() * springJocobian * leftVertexU;

      reducedSystemMatrix.block(colStart, colStart, colRank, colRank) += rightVertexU.transpose() * springJocobian * rightVertexU;

      reducedSystemMatrix.block(rowStart, colStart, rowRank, colRank) -= leftVertexU.transpose() * springJocobian * rightVertexU;

      reducedSystemMatrix.block(colStart, rowStart, colRank, rowRank) -= rightVertexU.transpose() * springJocobian * leftVertexU;
      
    }
  }
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::removeZeroBlocks(const MATRIX& input, MATRIX& output)
{
  int outputSize = 0;
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_reducedDofs[x] > 0)
      outputSize += _tetMesh->partitionRank(x);
  }

  output.conservativeResize(outputSize, outputSize);
  output.setZero();

  int outputRowStart = 0;
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_reducedDofs[x] == 0)
      continue;

    int outputColStart = 0;
    for(int y = 0; y < _tetMesh->totalPartitions(); y++){
      if(_reducedDofs[y] == 0)
        continue;
      output.block(outputRowStart, outputColStart, _tetMesh->partitionRank(x), _tetMesh->partitionRank(y)) = input.block(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRankStartIdx(y), _tetMesh->partitionRank(x), _tetMesh->partitionRank(y));
      
      outputColStart += _tetMesh->partitionRank(y);
    }
    outputRowStart += _tetMesh->partitionRank(x);
  }
}
template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::fillInZeroVectors(const VECTOR& input, VECTOR& output)
{
  output.conservativeResize(_tetMesh->totalPartitionRank());
  int inputOffset = 0;
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_reducedDofs[x] == 0){
      output.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)).setZero();
    }else{
      output.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) = input.segment(inputOffset, _tetMesh->partitionRank(x));
      inputOffset += _tetMesh->partitionRank(x);
    }
  }
}

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
void PARTITIONED_HYBRID_INTEGRATOR<FULL_MATERIAL_CACHE, SUB_MATERIAL_CACHE, BONE>::removeZeroVectors(const VECTOR& input, VECTOR& output)
{
  int outputSize = 0;
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_reducedDofs[x] > 0)
      outputSize += _tetMesh->partitionRank(x);
  }

  output.conservativeResize(outputSize);
  output.setZero();

  int outputOffset = 0;
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_reducedDofs[x] == 0)
      continue;
    output.segment(outputOffset, _tetMesh->partitionRank(x)) = input.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));
    outputOffset += _tetMesh->partitionRank(x);
  }
}
