#ifndef PARTITIONED_SUBSPACE_INTEGRATOR_H
#define PARTITIONED_SUBSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <iostream>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/SELF_COLLISION_DETECTOR.h>
#include <linearalgebra/LU_SOLVER.h>

using namespace::std;

template<class SUBSPACE_MATERIAL_CACHE, class BONE>
class PARTITIONED_SUBSPACE_INTEGRATOR
{
public:
  PARTITIONED_SUBSPACE_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~PARTITIONED_SUBSPACE_INTEGRATOR();

  const Real& energy()  { return _energy; };
  VECTOR& gradient()    { return _gradient; };
  MATRIX& hessian() { return _hessian; };
  LU_SOLVER& hessianInv() { return _hessianInv; };
  VECTOR& getPosition() { return _tetMesh->q(); };

  void initializeImplicitStep();

  void finalizeImplicitStep();

  MATRIX& computeSystemMatrix();
  
  Real computeSystemEnergy();
  
  VECTOR& computeSystemForce();
  
  void computeCollisionMatrices();
  
  void computeMaterialCache();

  SELF_COLLISION_DETECTOR<BONE>* scd() { return _scd; };
  Real computeSelfCollisionSpringForces(VECTOR& forceVector);
  Real computeExternalCollisionForces(VECTOR& forceVector);

  void computeSelfCollisionSpringForceJacobian(MATRIX& systemMatrix);
  void computeExternalCollisionForceJacobian(MATRIX& systemMatrix);

  void setPosition(VECTOR& newPosition);

  void resetState()
  {
    _tetMesh->q() = _initQ;

    _tetMesh->x() = _initX;

    _tetMesh->updateFullMesh();

    if(_dynamic)
      _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }

private:
  void computeInterfaceSpringJacobian(MATRIX& systemMatrix);
  Real computeInterfaceSpringForces(VECTOR& forceVector);

private:
  SUBSPACE_TET_MESH* _tetMesh;
  SUBSPACE_MATERIAL_CACHE* _subMaterialCache;

  SELF_COLLISION_DETECTOR<BONE>* _scd;
  RIGGER<BONE>* _rigger;

  VECTOR _tmpWorkspace;
  vector<VECTOR> _tmpWorkspace2;

  VECTOR _initX;
  VECTOR _initQ;

  VECTOR _gradient;

  MATRIX _stiffness;
  MATRIX _collisionJacobian;

  MATRIX _hessian;

  LU_SOLVER _hessianInv;

  // vector<MATRIX3> _transform;
  // vector<MATRIX3> _inverseTransform;

  vector<VECTOR> _diagMasses;

  // VECTOR _skinningDisp;
  Real _energy;

  bool _dynamic;

  MATRIX _reducedMass;

  VECTOR _positionOld;
  VECTOR _acceleration;
  VECTOR _velocityOld;

  MATRIX _fixedHessian;

  // dynamic integration timestep
  Real _dt;

  // openmp support
  MATRIX* _SCJacobianCopies;

private:
  Real _scfMultiplier;
  Real _interfaceSpringConstant;
};

#include "PARTITIONED_SUBSPACE_INTEGRATOR.inl"
#endif
