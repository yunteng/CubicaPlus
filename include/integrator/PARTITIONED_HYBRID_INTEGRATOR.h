#ifndef PARTITIONED_HYBRID_INTEGRATOR_H
#define PARTITIONED_HYBRID_INTEGRATOR_H

#include <SETTINGS.h>
#include <iostream>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/SELF_COLLISION_DETECTOR.h>
#include <linearalgebra/LU_SOLVER.h>

using namespace::std;

template<class FULL_MATERIAL_CACHE, class SUB_MATERIAL_CACHE, class BONE>
class PARTITIONED_HYBRID_INTEGRATOR
{
public:
  PARTITIONED_HYBRID_INTEGRATOR(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~PARTITIONED_HYBRID_INTEGRATOR();

  const Real& energy()  { return _energy; };
  VECTOR& gradient()    { return _gradient; };
  SpMat& hessian() { return _hessian; };

  void computeMaterialCache();

  void initializeImplicitStep();

  void finalizeImplicitStep();
  
  Real computeSystemEnergy();

  VECTOR& computeSystemForce();

  BLOCK_SPARSE_MATRIX& computeSystemMatrix();

  SELF_COLLISION_DETECTOR<BONE>*& scd() { return _scd; };
  Real computeSelfCollisionSpringForces(VECTOR& fullForceVector, VECTOR& reducedForceVector);
  Real computeExternalCollisionForces(VECTOR& forceVector);
  
  void computeSelfCollisionSpringForceJacobian(BLOCK_COO_MATRIX& fullSystemMatrix, MATRIX& reducedSystemMatrix);

  void computeCollisionMatrices();
  
  void computeExternalCollisionForceJacobian(BLOCK_COO_MATRIX& systemMatrix);

  void setPosition(VECTOR& newPosition);

  void computeQi(const VECTOR& u_s);
  
  void getMatVecMult(const VECTOR& input, VECTOR& output);

  int activeFullDofs() { return _tetMesh->partitionFullsimDofStartIdx(_tetMesh->totalPartitions()); };

  VECTOR& hessianDiagonal() { return _hessianDiagonal; };

  LU_SOLVER& reducedKiiInv() { return _reducedKiiInv; };


  void resetState()
  {
    _tetMesh->x() = _initPosition;
    _tetMesh->updateFullMesh();
    _tetMesh->changeToPartitionOrder(_initPosition, _partitionedX);
    if(_dynamic)
      _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }

  inline bool hasEntireFullPartitions() { return _hasEntireFullPartitions; }

  void removeZeroBlocks(const MATRIX& input, MATRIX& output);
  void fillInZeroVectors(const VECTOR& input, VECTOR& output);
  void removeZeroVectors(const VECTOR& input, VECTOR& output);

  VECTOR& isFullsim() { return _isFullsim; };

private:

  Real computeInterfaceSpringForces(VECTOR& fullForceVector, VECTOR& reducedForceVector);

  void computeInterfaceSpringJacobians(BLOCK_COO_MATRIX& fullSystemMatrix, MATRIX& reducedSystemMatrix);

  void computeReducedInterfaceSpringJacobians(MATRIX& reducedSystemMatrix);

  void computeFixedSystemMatrices();

  void computeFullsimRegions();

private:
  SUBSPACE_TET_MESH* _tetMesh;
  FULL_MATERIAL_CACHE* _fullMaterialCache;
  SUB_MATERIAL_CACHE* _subMaterialCache;
  SELF_COLLISION_DETECTOR<BONE>* _scd;
  RIGGER<BONE>* _rigger;

  vector<int> _previousFullDofs;
  vector<int>& _fullDofs;
  vector<int>& _reducedDofs;
  bool _hasEntireFullPartitions;

  VECTOR _qi;

  VECTOR _gradient;
  
  VECTOR _fullRegionGradient;
  VECTOR _reducedRegionGradient;

  SpMat _hessian;
  VECTOR _partitionedX;

  MATRIX _reducedKii;

  // in case there are zero diagonal blocks in _reducedKii, remove them
  MATRIX _prunedReducedKii;

  LU_SOLVER _reducedKiiInv;

  MATRIX _completeReducedInterfaceJacobians;

  vector<MATRIX> _transformedUs;
  vector<MATRIX> _untransformedUs;

  // vector<MATRIX3> _transform;
  // vector<MATRIX3> _inverseTransform;
  // VECTOR _skinningDisp;
  VECTOR _restDisp;
  VECTOR _workspace;
  VECTOR _workspace2;

  Real _energy;

  bool _dynamic;

  VECTOR _natualOrderedMass;
  VECTOR _diagMass;
  vector<MATRIX> _reducedMasses;

  BLOCK_COO_MATRIX _SCJacobians;
  MATRIX _reducedSCJacobians;

  SpMat _fixedFullHessian;
  SpMat _SCJacobiansSpMat;
  VECTOR _hessianDiagonal;
  VECTOR _stiffnessDiagonal;

  SpMat _stiffness;

  MATRIX _fixedReducedHessian;

  vector<SpMat> _UTKis;

  VECTOR _positionOld;
  VECTOR _velocity;
  VECTOR _acceleration;
  VECTOR _velocityOld;
  vector<VECTOR> _accelerationForces;

  VECTOR _initPosition;

  VECTOR _partitionedVelocity;
  VECTOR _keepDynamicRegion;

  Real _dt;

  VECTOR _isFullsim;

  vector<int> _previousCollisionPoints;

  // openmp support
  MATRIX* _SCJacobianCopies;

private:
  Real _scfMultiplier;
  Real _interfaceSpringConstant;
};

#include "PARTITIONED_HYBRID_INTEGRATOR.inl"
#endif
