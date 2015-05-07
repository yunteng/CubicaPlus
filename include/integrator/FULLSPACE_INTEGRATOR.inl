#include <util/SIMPLE_PARSER.h>
#include <util/MATRIX_UTIL.h>

template<class MATERIAL_CACHE, class BONE>
FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::FULLSPACE_INTEGRATOR(TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger),
  _scd(NULL)
{
  _materialCache = new MATERIAL_CACHE(tetMesh);
  _scd = new SELF_COLLISION_DETECTOR<BONE>(tetMesh, rigger);
  _scfMultiplier = SIMPLE_PARSER::getFloat("scf multiplier", 1.0);
  if(SIMPLE_PARSER::getBool("verbose", true))
    cout << " SCF multiplier: " << _scfMultiplier << endl;

  _dynamic = SIMPLE_PARSER::getBool("dynamic", false);
  if(_dynamic){
    _dt = SIMPLE_PARSER::getFloat("timestep", 1.0 / 60.0);
    if(SIMPLE_PARSER::getBool("verbose", true)){
      cout << " Dynamic simulation " << endl
           << "   timestep:      " << _dt << endl;
    }

    _tetMesh->recoverX();
    _positionOld = _tetMesh->x();
    _velocity.resize(_tetMesh->dofs());
    _acceleration.resize(_tetMesh->dofs());
    _velocityOld.resize(_tetMesh->dofs());
    
    _velocity.setZero();
    _acceleration.setZero();
    _velocityOld.setZero();
  }
}

template<class MATERIAL_CACHE, class BONE>
FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::~FULLSPACE_INTEGRATOR()
{
  if(_materialCache){
    delete _materialCache;
    _materialCache = NULL;
  }
  if(_scd){
    delete _scd;
    _scd = NULL;
  }
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::initializeImplicitStep()
{
  _tetMesh->recoverX();
  _initPosition = _tetMesh->x();
  
  _scd->vertexVsTetSCD();

  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();
  _individualCollisionResponses.resize(selfCollisionPoints.size());

  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){
    SELF_COLLISION_INFO& info = selfCollisionPoints[x];
    _individualCollisionResponses[x].vertexID = info.vertexID;
    _individualCollisionResponses[x].surfaceID = info.faceID;
    _individualCollisionResponses[x].vertexPartition = info.vertexPartition;
    _individualCollisionResponses[x].trianglePartition = info.trianglePartition;
    _individualCollisionResponses[x].collisionResponses.clear();
  }

  if(_dynamic)
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::finalizeImplicitStep()
{
  _tetMesh->updateFullMesh();

  if(_dynamic){
    _velocityOld = (1.0 / _dt) * (_tetMesh->x() - _positionOld);
    _positionOld = _tetMesh->x();
  }
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::setPosition(VECTOR& newPosition)
{
  _tetMesh->x() = newPosition;
  if(_dynamic){
    _acceleration = (1.0 / _dt / _dt) * (_tetMesh->x() - _positionOld) - (1.0 / _dt) * _velocityOld;
  }
}

template<class MATERIAL_CACHE, class BONE>
Real FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemEnergy()
{
  _energy = _materialCache->computeElasticEnergy();

  if(_dynamic)
    _energy += 0.5 * _dt * _dt * _acceleration.dot(_tetMesh->massMatrix() * _acceleration);

  return _energy;
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeCollisionMatrices()
{
  _collisionJacobian.conservativeResize(_hessian.rows(), _hessian.cols());
  _collisionJacobian.clear();

  computeSelfCollisionSpringForceJacobian(_collisionJacobian);
  computeExternalCollisionForceJacobian(_collisionJacobian);

  _collisionJacobian.order();
  _collisionJacobian.aggregate();

  if(_collisionJacobian.nnZ() > 0)
    _systemMatrixDiag = _hessian.diag() + _collisionJacobian.diag();
}

template<class MATERIAL_CACHE, class BONE>
COO_MATRIX& FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemMatrix()
{
  COO_MATRIX& stiffness = _materialCache->computeStiffnessMatrix();
  _hessian = stiffness;

  if(_dynamic){
    VECTOR& MVec = _tetMesh->massVector();
    for(int x = 0; x < MVec.size(); x++){
      Real val = (1.0 / _dt / _dt) * MVec[x];
      _hessian.add(val, x * 3, x * 3);
      _hessian.add(val, x * 3 + 1, x * 3 + 1);
      _hessian.add(val, x * 3 + 2, x * 3 + 2);
    }
  }

  computeCollisionMatrices();

  _hessian.order();
  _hessian.aggregate(false);

  _systemMatrixDiag = _hessian.diag();

  return _hessian;
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::
  computeMaterialCache()
{
  _materialCache->cacheDecompositions();
}

template<class MATERIAL_CACHE, class BONE>
VECTOR& FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSystemForce()
{
  VECTOR& R = _materialCache->computeInternalForce();

  _gradient = R * -1.0;

  _energy += computeSelfCollisionSpringForces(_gradient);

  _energy += computeExternalCollisionForces(_gradient);

  if(_dynamic)
    _gradient += _tetMesh->massMatrix() * _acceleration;

  return _gradient;
}

template<class MATERIAL_CACHE, class BONE>
Real FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForces(VECTOR& forceVector)
{
  Real energy = 0;

  vector<SELF_COLLISION_INFO>& selfCollisionPoints = _scd->selfCollisionPoints();

  for(unsigned int x = 0; x < selfCollisionPoints.size(); x++){

    VECTOR collisionResponses(12);

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

    VEC3F normal = surfacePosition - leftVertex;
    normal.normalize();

    MATRIX3 M = MATRIX_UTIL::outer_product(normal);
    M *= 0.5;
    M += 0.5 * MATRIX3::Identity();
    info.M = M;

    VEC3F penaltyForce = info.M * (surfacePosition - leftVertex) * multiplier;

    energy += 0.5 * (surfacePosition - leftVertex).dot(penaltyForce);

    if(!_tetMesh->isConstrained(leftIndex)){
      forceVector.segment<3>(leftIndex * 3) -= penaltyForce;
      collisionResponses.head<3>() = -penaltyForce;
    }
    for(int i = 0; i < 3; i++){
      if(!_tetMesh->isConstrained(triangleIndices[i])){
        forceVector.segment<3>(triangleIndices[i] * 3) += lambda[i] * penaltyForce;
        collisionResponses.segment<3>((i + 1) * 3) = lambda[i] * penaltyForce;
      }
    }
    _individualCollisionResponses[x].collisionResponses.push_back(collisionResponses);
  }

  return energy;
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeSelfCollisionSpringForceJacobian(COO_MATRIX& systemMatrix)
{
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

    VEC3F normal = surfacePosition - leftVertex;
    normal.normalize();

    MATRIX3 M = MATRIX_UTIL::outer_product(normal);
    M *= 0.5;
    M += 0.5 * MATRIX3::Identity();
    info.M = M;

    MATRIX3 leftJocobian = multiplier * M;

    if(!_tetMesh->isConstrained(leftIndex))
      systemMatrix.add3x3(leftJocobian, 3 * leftIndex, 3 * leftIndex);

    for(int x = 0; x < 3; x++){
      if(_tetMesh->isConstrained(triangleIndices[x]))
        continue;
      for(int y = 0; y < 3; y++){
        if(_tetMesh->isConstrained(triangleIndices[y]))
          continue;
        systemMatrix.add3x3(lambda[x] * lambda[y] * multiplier * M, triangleIndices[x] * 3, triangleIndices[y] * 3);
      }
    }
    if(_tetMesh->isConstrained(leftIndex))
      continue;

    for(int x = 0; x < 3; x++){
      Real J = -lambda[x] * multiplier;

      if(_tetMesh->isConstrained(triangleIndices[x]))
        continue;

      systemMatrix.add3x3(J * M, 3 * leftIndex, 3 * triangleIndices[x]);
      systemMatrix.add3x3(J * M.transpose(), 3 * triangleIndices[x], 3 * leftIndex);
    }
  }
}
template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeExternalCollisionForceJacobian(COO_MATRIX& systemMatrix)
{
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();

  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    MATRIX3 springJacobian = surface->springJacobian(*vertex);
    
    systemMatrix.add3x3(springJacobian, vertexID * 3, vertexID * 3);
  }
}
template<class MATERIAL_CACHE, class BONE>
Real FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::computeExternalCollisionForces(VECTOR& forceVector)
{
  Real energy = 0;
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs = _tetMesh->collisionPairs();
  for(unsigned int x = 0; x < collisionPairs.size(); x++){
    VEC3F* vertex = collisionPairs[x].first;
    SURFACE* surface = collisionPairs[x].second;
    int vertexID = _tetMesh->vertexID(vertex);

    VEC3F force = surface->force(*vertex);
    forceVector.segment<3>(vertexID * 3) += force;
    energy += 0.5 * force.squaredNorm() / surface->collisionStiffness();
  }
  return energy;
}

template<class MATERIAL_CACHE, class BONE>
void FULLSPACE_INTEGRATOR<MATERIAL_CACHE, BONE>::writeSelfCollisionResponses(const string& filename)
{
  int size = _individualCollisionResponses.size();
  if(size == 0)
    return;

  FILE* file = fopen(filename.c_str(), "wb");
  if(file == NULL){
    cout << __FILE__ << " " << __LINE__ << endl;
    cout << "cannot open " << filename << " to write!!!!" << endl;
    return;
  }
  fwrite((void*)&size, sizeof(int), 1, file);

  for(int x = 0; x < size; x++){
    int vertexPartition = _individualCollisionResponses[x].vertexPartition;
    int trianglePartition = _individualCollisionResponses[x].trianglePartition;
    int vertexID = _individualCollisionResponses[x].vertexID;
    int faceID = _individualCollisionResponses[x].surfaceID;
    
    vector<VECTOR>& collisionResponses = _individualCollisionResponses[x].collisionResponses;
    int n = collisionResponses.size();

    fwrite((void*)&vertexPartition, sizeof(int), 1, file);
    fwrite((void*)&trianglePartition, sizeof(int), 1, file);
    fwrite((void*)&vertexID, sizeof(int), 1, file);
    fwrite((void*)&faceID, sizeof(int), 1, file);
    
    fwrite((void*)&n, sizeof(int), 1, file);

    for(int y = 0; y < n; y++){
      for(int z = 0; z < 12; z++){
        Real val = collisionResponses[y][z];
        fwrite((void*)&val, sizeof(Real), 1, file);
      }
    }
  }
  fclose(file);
}
