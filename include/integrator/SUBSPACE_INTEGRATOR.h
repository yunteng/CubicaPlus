#ifndef SUBSPACE_INTEGRATOR_H
#define SUBSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <iostream>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/SELF_COLLISION_DETECTOR.h>
#include <linearalgebra/LU_SOLVER.h>

using namespace::std;

template<class MATERIAL_CACHE, class BONE>
class SUBSPACE_INTEGRATOR
{
public:
  SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~SUBSPACE_INTEGRATOR();

  const Real& energy()  { return _energy; };
  VECTOR& gradient()    { return _gradient; };
  MATRIX& hessian() { return _hessian; };
  LU_SOLVER& hessianInv() { return _hessianInv; };

  VECTOR& getPosition() { return _tetMesh->q(); };

  void initializeImplicitStep();

  void finalizeImplicitStep();
  
  Real computeSystemEnergy();

  VECTOR& computeSystemForce();

  MATRIX& computeSystemMatrix();

  void computeCollisionMatrices();

  inline MATRIX vertexTransformedU(int vertexID)
  {
    return _transformedU.block(vertexID * 3, 0, 3, _transformedU.cols());
  }
  
  void computeMaterialCache();

  SELF_COLLISION_DETECTOR<BONE>* scd() { return _scd; };
  Real computeSelfCollisionSpringForces(VECTOR& forceVector);

  void computeSelfCollisionSpringForceJacobian(MATRIX& systemMatrix);

  void setPosition(VECTOR& newPosition);

  void resetState()
  {
    _tetMesh->q() = _initPosition;

    if(_useTransformedBasis){
      _tetMesh->x() = _rigger->skinningDisp() + _transformedU * _tetMesh->q();
    }else{
      _tetMesh->x() = _tetMesh->U() * _tetMesh->q();
    }

    _tetMesh->updateFullMesh();
    if(_dynamic)
      _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }

private:
  SUBSPACE_TET_MESH* _tetMesh;
  MATERIAL_CACHE* _materialCache;

  SELF_COLLISION_DETECTOR<BONE>* _scd;
  RIGGER<BONE>* _rigger;

  VECTOR _gradient;
  MATRIX _hessian;

  LU_SOLVER _hessianInv;

  MATRIX _stiffness;
  MATRIX _collisionJacobian;
  
  VECTOR _diagMass;

  VECTOR _initPosition;

  Real _energy;

  bool _useKrysl;
  bool _useTransformedBasis;

  bool _dynamic;

  MATRIX _reducedMass;
  MATRIX _transformedU;
  VECTOR _positionOld;
  VECTOR _acceleration;
  VECTOR _velocityOld;

  // dynamic integration timestep
  Real _dt;

  // openmp support
  MATRIX* _SCJacobianCopies;

private:
  Real _scfMultiplier;
};

#include "SUBSPACE_INTEGRATOR.inl"
#endif
