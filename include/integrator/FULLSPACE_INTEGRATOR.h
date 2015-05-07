#ifndef FULLSPACE_INTEGRATOR_H
#define FULLSPACE_INTEGRATOR_H

#include <SETTINGS.h>
#include <iostream>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/TET_MESH.h>
#include <geometry/SELF_COLLISION_DETECTOR.h>

using namespace::std;

template<class MATERIAL_CACHE, class BONE>
class FULLSPACE_INTEGRATOR
{
public:
  FULLSPACE_INTEGRATOR(TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~FULLSPACE_INTEGRATOR();

  const Real& energy()  { return _energy; };
  VECTOR& gradient()    { return _gradient; };
  COO_MATRIX& hessian() { return _hessian; };
  VECTOR& getPosition() { return _tetMesh->x(); };

  void initializeImplicitStep();

  void finalizeImplicitStep();

  Real computeSystemEnergy();

  VECTOR& computeSystemForce();

  COO_MATRIX& computeSystemMatrix();
  
  void computeMaterialCache();

  void computeCollisionMatrices();

  SELF_COLLISION_DETECTOR<BONE>*& scd() { return _scd; };
  Real computeSelfCollisionSpringForces(VECTOR& forceVector);
  void computeSelfCollisionSpringForceJacobian(COO_MATRIX& systemMatrix);

  void setPosition(VECTOR& newPosition);

  void getMatVecMult(const VECTOR& input, VECTOR& output)
  { 
    output = _hessian * input;
    if(_collisionJacobian.nnZ() > 0)
      output += _collisionJacobian * input; 
  }
  
  void computeExternalCollisionForceJacobian(COO_MATRIX& systemMatrix);
  Real computeExternalCollisionForces(VECTOR& forceVector);
  void writeSelfCollisionResponses(const string& filename);

  static void staticGetMatVecMult(void* integrator, const VECTOR& input, VECTOR& output)
  {
    static_cast<FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>*>(integrator)->getMatVecMult(input, output);
  }
  VECTOR& systemMatrixDiag() { return _systemMatrixDiag; };

  void resetState()
  {
    _tetMesh->x() = _initPosition;
    _tetMesh->updateFullMesh();
    if(_dynamic)
      _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }
  
private:
  TET_MESH* _tetMesh;
  MATERIAL_CACHE* _materialCache;
  SELF_COLLISION_DETECTOR<BONE>* _scd;
  RIGGER<BONE>* _rigger;

  VECTOR _gradient;
  COO_MATRIX _hessian;
  COO_MATRIX _collisionJacobian;
  Real _energy;

  VECTOR _systemMatrixDiag;
  VECTOR _initPosition;

  bool _dynamic;

  VECTOR _positionOld;
  VECTOR _velocity;
  VECTOR _acceleration;
  VECTOR _velocityOld;

  Real _dt;

  struct SELF_COLLISION_PAIR{
    int vertexID;
    int surfaceID;
    int vertexPartition;
    int trianglePartition;
    vector<VECTOR> collisionResponses;
  };
  vector<SELF_COLLISION_PAIR> _individualCollisionResponses;

private:
  Real _scfMultiplier;
};

#include "FULLSPACE_INTEGRATOR.inl"
#endif