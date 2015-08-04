/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
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

  LU_SOLVER& reducedKssInv() { return _reducedKssInv; };


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

  /*
  the dofs in each partition that will be 
  simulated in full space
  */
  vector<int>& _fullDofs;
  /*
  fullDofs in the previous frame
  */
  vector<int> _previousFullDofs;
  /*
  the dofs in each partition that will be
  simulated in subspace
  */
  vector<int>& _reducedDofs;

  /*
  if any partition is entirely simulated 
  in full space
  */
  bool _hasEntireFullPartitions;

  /*
  the delta subspace coordinates
  */
  VECTOR _qi;

  VECTOR _gradient;
  
  VECTOR _fullRegionGradient;
  VECTOR _reducedRegionGradient;

  SpMat _hessian;

  /*
  partitioned world space displacement vector
  */
  VECTOR _partitionedX;

  MATRIX _reducedKss;

  // in case there are zero diagonal blocks in _reducedKss, remove them
  MATRIX _prunedReducedKss;

  LU_SOLVER _reducedKssInv;

  MATRIX _reducedInterfaceJacobians;

  /*
  transformed (with skinning rotation) basis
  the rows are reordered according to current
  condensation partitioning
  */
  vector<MATRIX> _transformedUs;
  /*
  untransformed basis, the rows are reordered 
  according to current condensation partitioning
  */
  vector<MATRIX> _untransformedUs;

  VECTOR _workspace;
  VECTOR _workspace2;

  Real _energy;

  bool _dynamic;

  VECTOR _defaultOrderMass;
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

  /*
  U_s^T * K_sf
  */
  vector<SpMat> _UTKsf;

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
  // self collision penalty spring constant
  Real _scfMultiplier;
  // partition interface spring constant
  Real _interfaceSpringConstant;
};

#include "PARTITIONED_HYBRID_INTEGRATOR.inl"
#endif
