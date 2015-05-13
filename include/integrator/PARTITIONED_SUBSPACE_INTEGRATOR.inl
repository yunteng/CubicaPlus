#include <util/SIMPLE_PARSER.h>
#include <util/TIMING_BREAKDOWN.h>
#include <util/MATRIX_UTIL.h>
#if USING_SUBSPACE_OPENMP
#include <omp.h>
#endif

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::PARTITIONED_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger)
{
  _subMaterialCache = new SUBSPACE_MATERIAL_CACHE(tetMesh);

  _scd = new SELF_COLLISION_DETECTOR<BONE>(tetMesh, rigger);
  _scfMultiplier = SIMPLE_PARSER::getFloat("scf multiplier", 1.0);
  _interfaceSpringConstant = SIMPLE_PARSER::getFloat("interface spring constant", 100.0);

  if(SIMPLE_PARSER::getBool("verbose", true)){
    cout << " SCF multiplier: " << _scfMultiplier << endl;
    cout << " Interface spring constant: " << _interfaceSpringConstant << endl;
  }

  _SCJacobianCopies = new MATRIX[_tetMesh->totalCores()];
  for(int x = 0; x < _tetMesh->totalCores(); x++){
    _SCJacobianCopies[x].resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  }

  _fixedHessian.resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  _fixedHessian.setZero();
  computeInterfaceSpringJacobian(_fixedHessian);

  _dynamic = SIMPLE_PARSER::getBool("dynamic", false);

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

    VECTOR natualOrderedMass(_tetMesh->dofs());
    VECTOR& MVec = _tetMesh->massVector();
    for(int x = 0; x < MVec.size(); x++){
      natualOrderedMass[x * 3] = MVec[x];
      natualOrderedMass[x * 3 + 1] = MVec[x];
      natualOrderedMass[x * 3 + 2] = MVec[x];
    }

    _tetMesh->changeToPartitionOrder(natualOrderedMass, _diagMasses);

    _reducedMass.resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
    _reducedMass.setZero();
    for(unsigned int x = 0; x < _diagMasses.size(); x++){
      _reducedMass.block(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x), _tetMesh->partitionRank(x)) = _tetMesh->partitionBasis(x).transpose() * _diagMasses[x].asDiagonal() * _tetMesh->partitionBasis(x);
    }
    _fixedHessian += (1.0 / _dt / _dt) * _reducedMass;
    // TODO::start the simulation from non rest pose
  }
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::~PARTITIONED_SUBSPACE_INTEGRATOR()
{
  delete[] _SCJacobianCopies;

  if(_subMaterialCache)
    delete _subMaterialCache;
  if(_scd)
    delete _scd;
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::initializeImplicitStep()
{
  TIMING_BREAKDOWN::tic();
  _tetMesh->recoverX();
  TIMING_BREAKDOWN::toc("Recover X");

  TIMING_BREAKDOWN::tic();
  _scd->vertexVsTetSCD();
  TIMING_BREAKDOWN::toc("Self Collision Detection");

  // recover q

  TIMING_BREAKDOWN::tic();
  VECTOR restDisp = _tetMesh->x() - _rigger->skinningDisp();
  vector<MATRIX3>& skinningRotation = _rigger->skinningRotation();

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
    _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) = _tetMesh->partitionBasis(x).transpose() * perPartitionRestDisp[x];
  }

  _initX = _tetMesh->x();
  _initQ = _tetMesh->q();

  TIMING_BREAKDOWN::toc("Recover q");

  TIMING_BREAKDOWN::tic();
  _subMaterialCache->cacheKeyTetTransforms(skinningRotation);
  TIMING_BREAKDOWN::toc("Cache Key Tet Transforms");
  
  if(_dynamic)
  {
    TIMING_BREAKDOWN::tic();
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;

    TIMING_BREAKDOWN::toc("compute acceleration");
  }

}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::finalizeImplicitStep()
{
  _tetMesh->updateFullMesh();
  if(_dynamic){
    _velocityOld = (1.0 / _dt) * (_tetMesh->x() - _positionOld);
    _positionOld = _tetMesh->x();
  }
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::setPosition(VECTOR& newPosition)
{
  _tetMesh->q() = newPosition;

  vector<VECTOR> restDisp(_tetMesh->totalPartitions());

  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    restDisp[x] = _tetMesh->partitionBasis(x) * _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x));
  }

  _tetMesh->restoreDefaultOrder(restDisp, _tmpWorkspace);

  VECTOR& worldDisp = _tetMesh->x();
  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _tetMesh->dofs() / 3; x++)
  {
    worldDisp.segment<3>(x * 3) = _rigger->skinningRotation()[x] * _tmpWorkspace.segment<3>(x * 3);
  }
  worldDisp += _rigger->skinningDisp();

  if(_dynamic){
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
Real PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeSystemEnergy()
{
  _energy = _subMaterialCache->computeElasticEnergy();
  
  if(_dynamic)
    _energy += 0.5 * _dt * _dt * _acceleration.dot(_tetMesh->massMatrix() * _acceleration);

  return _energy;
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
MATRIX& PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeSystemMatrix()
{
  _stiffness = _subMaterialCache->computeReducedStiffnessMatrix();
  _stiffness += _fixedHessian;

  computeCollisionMatrices();

  return _hessian;
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeCollisionMatrices()
{
  _collisionJacobian.conservativeResize(_stiffness.rows(), _stiffness.cols());
  _collisionJacobian.setZero();

  computeSelfCollisionSpringForceJacobian(_collisionJacobian);
  computeExternalCollisionForceJacobian(_collisionJacobian);

  _hessian = _stiffness + _collisionJacobian;
  _hessianInv.compute(_hessian);
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::
  computeMaterialCache()
{
  _subMaterialCache->cacheDecompositions();
}
template<class SUBSPACE_MATERIAL_CACHE, class BONE>
VECTOR& PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeSystemForce()
{
  VECTOR& R = _subMaterialCache->computeInternalForce();

  _gradient = R * -1;
  _energy += computeSelfCollisionSpringForces(_gradient);
  _energy += computeExternalCollisionForces(_gradient);
  _energy += computeInterfaceSpringForces(_gradient);

  if(_dynamic){
    _tmpWorkspace.conservativeResize(_acceleration.size());

    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(unsigned int x = 0; x < _acceleration.size() / 3; x++)
    {
      _tmpWorkspace.segment<3>(x * 3) = _rigger->skinningRotation()[x].transpose() * _acceleration.segment<3>(x * 3);
    }
    _tmpWorkspace = _tetMesh->massMatrix() * _tmpWorkspace;
    _tmpWorkspace *= 1.0 / _dt / _dt;

    _tetMesh->changeToPartitionOrder(_tmpWorkspace, _tmpWorkspace2);

    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int x = 0; x < _tetMesh->totalPartitions(); x++)
    {
      _gradient.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) += _tetMesh->partitionBasis(x).transpose() * _tmpWorkspace2[x];
    }
  }

  return _gradient;
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
Real PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForces(VECTOR& forceVector)
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
    TODO::add dynamic
    */

    energy += 0.5 * (surfacePosition - leftVertex).dot(penaltyForce);

    if(!_tetMesh->isConstrained(leftIndex)){
      pair<int, int>& partitionID = _tetMesh->partitionedVertexID(leftIndex);
      forceVector.segment(_tetMesh->partitionRankStartIdx(partitionID.first), _tetMesh->partitionRank(partitionID.first))
       -= _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * (_rigger->skinningRotation()[leftIndex].transpose() * penaltyForce);
    }

    for(int i = 0; i < 3; i++){
      if(!_tetMesh->isConstrained(triangleIndices[i])){

        pair<int, int>& partitionID = _tetMesh->partitionedVertexID(triangleIndices[i]);

        forceVector.segment(_tetMesh->partitionRankStartIdx(partitionID.first), _tetMesh->partitionRank(partitionID.first))
         += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * (_rigger->skinningRotation()[triangleIndices[i]].transpose() * (lambda[i] * penaltyForce));
      }
    }
  }
  return energy;
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForceJacobian(MATRIX& systemMatrix)
{
  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();
  if(selfCollisionPoints.size() == 0)
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
    #pragma omp for schedule(static)
    #endif
    for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){

      SELF_COLLISION_INFO& info = selfCollisionPoints[x];

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

      /*
      TODO::add dynamic
      */

      if(!_tetMesh->isConstrained(leftIndex)){
        leftJocobian = _rigger->skinningRotation()[leftIndex].transpose() * leftJocobian * _rigger->skinningRotation()[leftIndex];

        pair<int, int>& partitionID = _tetMesh->partitionedVertexID(leftIndex);

        int rowStart = _tetMesh->partitionRankStartIdx(partitionID.first);
        int colStart = rowStart;
        int subRank = _tetMesh->partitionRank(partitionID.first);

        _SCJacobianCopies[id].block(rowStart, colStart, subRank, subRank) += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * leftJocobian * _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second);
      }

      for(int x = 0; x < 3; x++){
        if(_tetMesh->isConstrained(triangleIndices[x]))
          continue;
        pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(triangleIndices[x]);
        int rowStart = _tetMesh->partitionRankStartIdx(leftPartitionID.first);
        int rowRank  = _tetMesh->partitionRank(leftPartitionID.first);

        for(int y = 0; y < 3; y++){
          if(_tetMesh->isConstrained(triangleIndices[y]))
            continue;

          pair<int, int>& rightPartitionID = _tetMesh->partitionedVertexID(triangleIndices[y]);
          int colStart = _tetMesh->partitionRankStartIdx(rightPartitionID.first);
          int colRank  = _tetMesh->partitionRank(rightPartitionID.first);

          MATRIX3 rightJacobian = lambda[x] * lambda[y] * multiplier * M;

          rightJacobian = _rigger->skinningRotation()[triangleIndices[x]].transpose() * rightJacobian * _rigger->skinningRotation()[triangleIndices[y]];

          _SCJacobianCopies[id].block(rowStart, colStart, rowRank, colRank) += 
            _tetMesh->partitionVertexBasis(leftPartitionID.first, leftPartitionID.second).transpose() * rightJacobian * _tetMesh->partitionVertexBasis(rightPartitionID.first, rightPartitionID.second);
        }
      }
      if(_tetMesh->isConstrained(leftIndex))
        continue;

      pair<int, int>& leftPartitionID = _tetMesh->partitionedVertexID(leftIndex);
      int rowStart = _tetMesh->partitionRankStartIdx(leftPartitionID.first);
      int rowRank = _tetMesh->partitionRank(leftPartitionID.first);

      for(int x = 0; x < 3; x++){
        if(_tetMesh->isConstrained(triangleIndices[x]))
          continue;

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
    systemMatrix += _SCJacobianCopies[x];
}
template<class SUBSPACE_MATERIAL_CACHE, class BONE>
Real PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeInterfaceSpringForces(VECTOR& forceVector)
{
  Real energy = 0;
  
  map<pair<int, int>, vector<MATRIX> >& interfaceUs = _tetMesh->interfaceUs();

  for(map<pair<int, int>, vector<MATRIX> >::iterator iter = interfaceUs.begin(); iter != interfaceUs.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;
    vector<MATRIX>& springUs = iter->second;
    VECTOR leftq = _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(leftPartition), _tetMesh->partitionRank(leftPartition));
    VECTOR rightq = _tetMesh->q().segment(_tetMesh->partitionRankStartIdx(rightPartition), _tetMesh->partitionRank(rightPartition));

    VECTOR leftDisp = springUs[0] * leftq;
    VECTOR rightDisp = springUs[1] * rightq;
    VECTOR diff = leftDisp - rightDisp;
    energy += diff.squaredNorm() * 0.5 * _interfaceSpringConstant;

    leftq *= _interfaceSpringConstant;
    rightq *= _interfaceSpringConstant;

    forceVector.segment(_tetMesh->partitionRankStartIdx(leftPartition), _tetMesh->partitionRank(leftPartition)) += springUs[2] * leftq - springUs[3] * rightq;

    forceVector.segment(_tetMesh->partitionRankStartIdx(rightPartition), _tetMesh->partitionRank(rightPartition)) -= springUs[3].transpose() * leftq - springUs[4] * rightq;
  }
  cout << "interface spring energy " << energy << endl;
  return energy;
}
template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeInterfaceSpringJacobian(MATRIX& systemMatrix)
{
  map<pair<int, int>, vector<pair<int, int> > > interfaceSprings = _tetMesh->interfaceSprings();

  MATRIX3 springJocobian = _interfaceSpringConstant * MATRIX3::Identity();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = interfaceSprings.begin(); iter != interfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;
    vector<pair<int, int> >& springs = iter->second;
    for(unsigned int x = 0; x < springs.size(); x++){
      int leftVID = springs[x].first;
      int rightVID = springs[x].second;
      int vertexID = _tetMesh->partitionedVertices(leftPartition)[leftVID];
      assert(vertexID == _tetMesh->partitionedVertices(rightPartition)[rightVID]);

      int rowStart = _tetMesh->partitionRankStartIdx(leftPartition);
      int colStart = _tetMesh->partitionRankStartIdx(rightPartition);
      int rowRank = _tetMesh->partitionRank(leftPartition);
      int colRank = _tetMesh->partitionRank(rightPartition);

      MATRIX leftVertexU = _tetMesh->partitionVertexBasis(leftPartition, leftVID);
      MATRIX rightVertexU = _tetMesh->partitionVertexBasis(rightPartition, rightVID);

      systemMatrix.block(rowStart, rowStart, rowRank, rowRank) += leftVertexU.transpose() * springJocobian * leftVertexU;

      systemMatrix.block(colStart, colStart, colRank, colRank) += rightVertexU.transpose() * springJocobian * rightVertexU;

      systemMatrix.block(rowStart, colStart, rowRank, colRank) -= leftVertexU.transpose() * springJocobian * rightVertexU;

      systemMatrix.block(colStart, rowStart, colRank, rowRank) -= rightVertexU.transpose() * springJocobian * leftVertexU;
    }
  }
}
template<class SUBSPACE_MATERIAL_CACHE, class BONE>
void PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeExternalCollisionForceJacobian(MATRIX& systemMatrix)
{
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();

  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    pair<int, int>& partitionID = _tetMesh->partitionedVertexID(vertexID);

    int rowStart = _tetMesh->partitionRankStartIdx(partitionID.first);
    int colStart = rowStart;
    int subRank = _tetMesh->partitionRank(partitionID.first);

    MATRIX3 springJacobian = surface->springJacobian(*vertex);

    springJacobian = _rigger->skinningRotation()[vertexID].transpose() * springJacobian * _rigger->skinningRotation()[vertexID];
    
    systemMatrix.block(rowStart, colStart, subRank, subRank) += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * springJacobian * _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second);
  }
}

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
Real PARTITIONED_SUBSPACE_INTEGRATOR<SUBSPACE_MATERIAL_CACHE, BONE>::computeExternalCollisionForces(VECTOR& forceVector)
{
  Real energy = 0;

  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();

  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    pair<int, int>& partitionID = _tetMesh->partitionedVertexID(vertexID);

    int rowStart = _tetMesh->partitionRankStartIdx(partitionID.first);
    
    int subRank = _tetMesh->partitionRank(partitionID.first);

    VEC3F force = surface->force(*vertex);

    energy += 0.5 * force.squaredNorm() / surface->collisionStiffness();
    
    forceVector.segment(rowStart, subRank) += _tetMesh->partitionVertexBasis(partitionID.first, partitionID.second).transpose() * (_rigger->skinningRotation()[vertexID] * force);
  }
  return energy;
}
