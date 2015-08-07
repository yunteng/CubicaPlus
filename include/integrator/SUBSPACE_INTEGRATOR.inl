/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <util/SIMPLE_PARSER.h>
#include <util/TIMING_BREAKDOWN.h>
#include <util/MATRIX_UTIL.h>
#if USING_SUBSPACE_OPENMP
#include <omp.h>
#endif

template<class MATERIAL_CACHE, class BONE>
SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger)
{
  _materialCache = new MATERIAL_CACHE(tetMesh);
  _scd = new SELF_COLLISION_DETECTOR<BONE>(tetMesh, rigger);
  _scfMultiplier = SIMPLE_PARSER::getFloat("scf multiplier", 1.0);
  if(SIMPLE_PARSER::getBool("verbose", true))
    cout << " SCF multiplier: " << _scfMultiplier << endl;

  _useKrysl = SIMPLE_PARSER::getBool("use krysl", false);
  _useTransformedBasis = _tetMesh->isSkinningBasis();

  if(_useTransformedBasis)
    _transformedU.resizeLike(_tetMesh->U());

  _SCJacobianCopies = new MATRIX[_tetMesh->totalCores()];
  for(int x = 0; x < _tetMesh->totalCores(); x++){
    _SCJacobianCopies[x].resize(_tetMesh->rank(), _tetMesh->rank());
  }

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

    _tetMesh->q().resize(_tetMesh->rank());
    _tetMesh->q().setZero();

    _diagMass.resize(_tetMesh->dofs());
    _diagMass.setZero();
    VECTOR& MVec = _tetMesh->massVector();
    for(int x = 0; x < MVec.size(); x++){
      _diagMass[x * 3] = MVec[x];
      _diagMass[x * 3 + 1] = MVec[x];
      _diagMass[x * 3 + 2] = MVec[x];
    }

    _reducedMass = _tetMesh->U().transpose() * _tetMesh->massMatrix().rightMult(_tetMesh->U());

    _reducedMass *= 1.0 / _dt / _dt;
  }
}

template<class MATERIAL_CACHE, class BONE>
SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::~SUBSPACE_INTEGRATOR()
{
  delete[] _SCJacobianCopies;

  if(_materialCache)
    delete _materialCache;
  if(_scd)
    delete _scd;
}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::initializeImplicitStep()
{
  TIMING_BREAKDOWN::tic();
  _tetMesh->recoverX();
  TIMING_BREAKDOWN::toc("Recover X");

  TIMING_BREAKDOWN::tic();
  _scd->vertexVsTetSCD();
  TIMING_BREAKDOWN::toc("Self Collision Detection");

  // recover q
  if(_useTransformedBasis){

    TIMING_BREAKDOWN::tic();
    VECTOR restDisp = _tetMesh->x() - _rigger->skinningDisp();
    vector<MATRIX3>& skinningRotation = _rigger->skinningRotation();

    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(unsigned int x = 0; x < skinningRotation.size(); x++){
      restDisp.segment<3>(x * 3) = skinningRotation[x].transpose() * restDisp.segment<3>(x * 3);
    }

    _tetMesh->q() = _tetMesh->U().transpose() * restDisp;
    TIMING_BREAKDOWN::toc("Recover q");

    if(!_useKrysl){
      TIMING_BREAKDOWN::tic();
      _materialCache->cacheKeyTetTransforms(skinningRotation);
      TIMING_BREAKDOWN::toc("Cache Key Tet Transforms");
    }


    TIMING_BREAKDOWN::tic();
    int cols = _tetMesh->U().cols();
    assert(_tetMesh->U().rows() > 0);
    assert(cols > 0);
    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int x = 0; x < _tetMesh->U().rows() / 3; x++){
      _transformedU.block(x * 3, 0, 3, cols) = _rigger->skinningRotation()[x] * _tetMesh->U().block(x * 3, 0, 3, cols);
    }
    TIMING_BREAKDOWN::toc("Compute transformed U");

    if(_dynamic && _rigger->skinningMethod().compare("dual quaternion")){
      TIMING_BREAKDOWN::tic();

      _reducedMass = _transformedU.transpose() * _diagMass.asDiagonal() * _transformedU;
      _reducedMass *= 1.0 / _dt / _dt;

      TIMING_BREAKDOWN::toc("Compute reduced Mass");
    }

  }else{
    TIMING_BREAKDOWN::tic();
    _tetMesh->q() = _tetMesh->U().transpose() * _tetMesh->x();
    TIMING_BREAKDOWN::toc("Recover q");

  }

  _initPosition = _tetMesh->q();
  
  if(_dynamic)
  {
    TIMING_BREAKDOWN::tic();
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
    TIMING_BREAKDOWN::toc("compute acceleration");
  }

}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::finalizeImplicitStep()
{
  _tetMesh->updateFullMesh();
  if(_dynamic){
    _velocityOld = (1.0 / _dt) * (_tetMesh->x() - _positionOld);
    _positionOld = _tetMesh->x();
  }
}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::setPosition(VECTOR& newPosition){
  _tetMesh->q() = newPosition;

  if(_useTransformedBasis){
    _tetMesh->x() = _rigger->skinningDisp() + _transformedU * _tetMesh->q();
  }else{
    _tetMesh->x() = _tetMesh->U() * _tetMesh->q();
  }
  if(_dynamic){
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }
}

template<class MATERIAL_CACHE, class BONE>
Real SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemEnergy()
{
  _energy = _materialCache->computeElasticEnergy();
  
  if(_dynamic)
    _energy += 0.5 * _dt * _dt * _acceleration.dot(_tetMesh->massMatrix() * _acceleration);

  return _energy;
}

template<class MATERIAL_CACHE, class BONE>
MATRIX& SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemMatrix()
{
  if(_useKrysl){
    _materialCache->computeStiffnessMatrix();

    if(_useTransformedBasis){
      _stiffness = _materialCache->reduceStiffness(_transformedU);
    }else{
      _stiffness = _materialCache->reduceStiffness(_tetMesh->U());
    }

  }else{
    _stiffness = _materialCache->computeReducedStiffnessMatrix();
  }
  if(_dynamic)
    _stiffness += _reducedMass;

  computeCollisionMatrices();

  return _hessian;
}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeCollisionMatrices()
{
  _collisionJacobian.conservativeResize(_stiffness.rows(), _stiffness.cols());
  _collisionJacobian.setZero();
  computeSelfCollisionSpringForceJacobian(_collisionJacobian);
  _hessian = _stiffness + _collisionJacobian;

  _hessianInv.compute(_hessian);
}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::
  computeMaterialCache()
{
  _materialCache->cacheDecompositions();
}
template<class MATERIAL_CACHE, class BONE>
VECTOR& SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemForce()
{
  VECTOR& R = _materialCache->computeInternalForce();
  if(_useKrysl){
    R *= -1;
    if(_dynamic)
      R += _diagMass.asDiagonal() * _acceleration;

    if(_useTransformedBasis){
      _gradient = _transformedU.transpose() * R;
    }else{
      _gradient = _tetMesh->U().transpose() * R;
    }

  }else{
    _gradient = R * -1;
    if(_dynamic){
      if(_useTransformedBasis){
        _gradient += _transformedU.transpose() * (_diagMass.asDiagonal() * _acceleration);
      }else{
        _gradient += _tetMesh->U().transpose() * (_diagMass.asDiagonal() * _acceleration);
      }
    }
  }
  _energy += computeSelfCollisionSpringForces(_gradient);

  return _gradient;
}

template<class MATERIAL_CACHE, class BONE>
Real SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForces(VECTOR& forceVector)
{
  Real energy = 0;
  
  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();

  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){
    SELF_COLLISION_INFO& info = selfCollisionPoints[x];

    int leftIndex = info.vertexID;

    VEC3F& leftVertex = _tetMesh->vertices()[leftIndex];
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

    /*
    TODO::add dynamic
    */

    energy += 0.5 * (surfacePosition - leftVertex).dot(penaltyForce);

    if(!_tetMesh->isConstrained(leftIndex)){
      if(_useTransformedBasis){
        forceVector -= vertexTransformedU(leftIndex).transpose() * penaltyForce;
      }else{
        forceVector -= _tetMesh->vertexBasis(leftIndex).transpose() * penaltyForce;
      }
    }

    for(int i = 0; i < 3; i++){
      if(!_tetMesh->isConstrained(triangleIndices[i])){
        if(_useTransformedBasis){
          forceVector += vertexTransformedU(triangleIndices[i]).transpose() * (lambda[i] * penaltyForce);
        }else{
          forceVector += _tetMesh->vertexBasis(triangleIndices[i]).transpose() * (lambda[i] * penaltyForce);
        }
      }
    }
  }

  return energy;
}

template<class MATERIAL_CACHE, class BONE>
void SUBSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForceJacobian(MATRIX& systemMatrix)
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

      /*
      TODO::add dynamic
      */

      if(!_tetMesh->isConstrained(leftIndex)){
        if(_useTransformedBasis){
          leftJocobian = _rigger->skinningRotation()[leftIndex].transpose() * leftJocobian * _rigger->skinningRotation()[leftIndex];
        }

        _SCJacobianCopies[id] += _tetMesh->vertexBasis(leftIndex).transpose() * leftJocobian * _tetMesh->vertexBasis(leftIndex);
      }

      for(int x = 0; x < 3; x++){
        if(_tetMesh->isConstrained(triangleIndices[x]))
          continue;
        for(int y = 0; y < 3; y++){
          if(_tetMesh->isConstrained(triangleIndices[y]))
            continue;

          MATRIX3 rightJacobian = lambda[x] * lambda[y] * multiplier * M;
          if(_useTransformedBasis){
            rightJacobian = _rigger->skinningRotation()[triangleIndices[x]].transpose() * rightJacobian * _rigger->skinningRotation()[triangleIndices[y]];
          }

          _SCJacobianCopies[id] += _tetMesh->vertexBasis(triangleIndices[x]).transpose() * rightJacobian * _tetMesh->vertexBasis(triangleIndices[y]);
        }
      }
      if(_tetMesh->isConstrained(leftIndex))
        continue;

      for(int x = 0; x < 3; x++){
        MATRIX jacobian = -lambda[x] * multiplier * M;
        if(_tetMesh->isConstrained(triangleIndices[x]))
          continue;
        if(_useTransformedBasis){
          jacobian = _rigger->skinningRotation()[leftIndex].transpose() * jacobian * _rigger->skinningRotation()[triangleIndices[x]];
        }

        _SCJacobianCopies[id] += _tetMesh->vertexBasis(leftIndex).transpose() * jacobian * _tetMesh->vertexBasis(triangleIndices[x]);

        _SCJacobianCopies[id] += _tetMesh->vertexBasis(triangleIndices[x]).transpose() * jacobian.transpose() * _tetMesh->vertexBasis(leftIndex);
      }
    }
  }
  for(int x = 0; x < _tetMesh->totalCores(); x++)
    systemMatrix += _SCJacobianCopies[x];
}
