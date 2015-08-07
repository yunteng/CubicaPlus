/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <geometry/TET_MESH.h>
#include <set>
#include <queue>

void TET_MESH::generateF()
{
  if(_F.size() != _tets.size())
    _F.resize(_tets.size());

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int x = 0; x < _tets.size(); x++){
    _F[x] = _tets[x].F();
  }
}

void TET_MESH::recoverX()
{
  assert(_unconstrainedSize * 3 == _x.size());
#ifdef USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _unconstrainedSize; i++)
  {
    _x.segment<3>(3 * i) = _vertices[i] - _restPose[i];
  }
  // cout << "done" << endl;
}

void TET_MESH::updateFullMesh()
{
  // cout << "update fullmesh" << endl;
  assert(_unconstrainedSize * 3 == _x.size());
#ifdef USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int i = 0; i < _unconstrainedSize; i++)
  {
    _vertices[i] = _restPose[i] + _x.segment<3>(3 * i);
  }

  // for(unsigned int x = 0; x < _clonedVertices.size(); x++){
  //   int vid = _clonedVertexOID[x];
  //   _clonedVertices[x] = _restPose[vid] + _clonedX.segment<3>(3 * x);
  // }
  // cout << "done" << endl;
}

void TET_MESH::updateFullMesh(const VECTOR& pos)
{
  assert(_unconstrainedSize * 3 == pos.size());
#ifdef USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _unconstrainedSize; i++)
  {
    _vertices[i] = _restPose[i] + pos.segment<3>(3 * i);
  }
}

void TET_MESH::reset()
{
  _x.setZero();
  updateFullMesh();
  for(unsigned int x = _unconstrainedSize; x < _vertices.size(); x++)
    _vertices[x] = _restPose[x];
}

void TET_MESH::smoothSurface()
{
  if(_surfaceVertexOneRings.empty())
    computeSurfaceVertexOneRings();

  vector<VEC3F> center(_surfaceVertices.size());
  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(unsigned int x = 0; x < _surfaceVertexOneRings.size(); x++)
  {
    center[x].setZero();
    Real sumWeight = 0;
    for(unsigned int y = 0; y < _surfaceVertexOneRings[x].size(); y++){
      Real dist = (*_surfaceVertices[x] - (*_surfaceVertexOneRings[x][y])).norm();
      sumWeight += dist;
      center[x] += *_surfaceVertexOneRings[x][y] * dist;
    }
    center[x] /= sumWeight;
  }
  Real strength = 1.0;
  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++){
    *_surfaceVertices[x] = strength * center[x] + (1 - strength) * (*_surfaceVertices[x]);
  }

}
void TET_MESH::computeSurfaceVertexNormal()
{
  _surfaceNormals.conservativeResize(_surfaceVertexSize * 3);
  _surfaceNormals.setZero();
  for(unsigned int x = 0; x < _surfaceFaces.size(); x++){
    TRIANGLE& face = _surfaceFaces[x];
    for(int j = 0; j < 3; j++){
      const VEC3F& v0 = *(face.vertices[(j + 3 - 1) % 3]);
      const VEC3F& v1 = *(face.vertices[j]);
      const VEC3F& v2 = *(face.vertices[(j + 1) % 3]);
      VEC3F e1 = (v0 - v1).normalized();
      VEC3F e2 = (v2 - v1).normalized();
      VEC3F normal = e1.cross(e2).normalized();
      Real cosAngle = e1.dot(e2);
      Real angle = acos(cosAngle);
      _surfaceNormals.segment<3>(_vertexID[face.vertices[j]] * 3) -= angle * normal;
    }
  }
  for(int x = 0; x < _surfaceVertexSize; x++){
    _surfaceNormals.segment<3>(x * 3) /= _surfaceNormals.segment<3>(x * 3).norm();
  }
}
void TET_MESH::computeSurfaceVertexOneRings()
{
  cout << "compute surface vertex one ring...";
  flush(cout);

  vector<set<int> > tempAdj(_vertices.size());
  for(unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE& face = _surfaceFaces[x];
    int indices[3];
    for(int y = 0; y < 3; y++)
      indices[y] = _vertexID[face.vertices[y]];
    for(int y = 0; y < 3; y++)
    {
      tempAdj[indices[y]].insert(indices[(y + 1) % 3]);
      tempAdj[indices[y]].insert(indices[(y + 2) % 3]);
    } 
  }
  _surfaceVertexOneRings.clear();
  _surfaceVertexOneRings.resize(_surfaceVertices.size());
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++)
  {
    int vertexID = _vertexID[_surfaceVertices[x]];
    _surfaceVertexOneRings[x].clear();
    for(set<int>::iterator iter = tempAdj[vertexID].begin(); iter != tempAdj[vertexID].end(); iter++)
      _surfaceVertexOneRings[x].push_back(&_vertices[*iter]);
  }
  cout << "done." << endl;
}

void TET_MESH::computeVertexAdjacency(vector<vector<int> >& adjacentVerts)
{
  vector<set<int> > tempAdj(_vertices.size());

  for(unsigned int x = 0; x < _tets.size(); x++){
    TET& tet = _tets[x];
    int indices[4];
    for(int y = 0; y < 4; y++)
      indices[y] = _vertexID[tet.vertices[y]];
    for(int y = 0; y < 4; y++){
      tempAdj[indices[y]].insert(indices[(y + 1) % 4]);
      tempAdj[indices[y]].insert(indices[(y + 2) % 4]);
      tempAdj[indices[y]].insert(indices[(y + 3) % 4]);
    }
  }
  adjacentVerts.resize(_vertices.size());
  for(unsigned int x = 0; x < tempAdj.size(); x++){
    
    copy(tempAdj[x].begin(), tempAdj[x].end(), back_inserter(adjacentVerts[x]));
  }
}

void TET_MESH::computeMasses()
{
  _masses.resize(3 * _unconstrainedSize, 3 * _unconstrainedSize);
  _massVec.resize(_unconstrainedSize);

  _totalVolume = 0;
  map<VEC3F*, Real> vertexMass;
  for(unsigned x = 0; x < _tets.size(); x++){
    Real volume = _tets[x].volume();
    for(int y = 0; y < 4; y++){
      vertexMass[_tets[x].vertices[y]] += volume / 4;
    }
    _totalVolume += volume;
  }
  Real massSum = 0;
  for(int x = 0; x < _unconstrainedSize; x++){
    Real myMass = vertexMass[&_vertices[x]] / _totalVolume;
    massSum += myMass;
    _massVec[x] = myMass;
    for(int y = 0; y < 3; y++)
      _masses.add(myMass, x * 3 + y, x * 3 + y);
  }
  Real totalMass = SIMPLE_PARSER::getFloat("total mass", 1.0);
  _masses *= totalMass;

  _massVec *= totalMass;

  _totalMass = totalMass * massSum;


  cout << " Total volume " << _totalVolume << endl;
  cout << " Sum of mass " << _totalMass << endl;
}

void TET_MESH::computeLaplacianMatrix(const vector<vector<int> >& adjacentVerts, SpMat& matrix)
{
  vector<VEC3F>& vertices = restPose();
  matrix.resize(vertices.size(), vertices.size());

  int maxRowEle = 0;
  for(unsigned int x = 0; x < adjacentVerts.size(); x++)
    if(adjacentVerts[x].size() > maxRowEle)
      maxRowEle = adjacentVerts[x].size();

  matrix.reserve(Eigen::VectorXi::Constant(vertices.size(), maxRowEle + 1));

  for(unsigned int x = 0; x < vertices.size(); x++){
    const VEC3F& vpos = vertices[x];
    Real sum = 0;
    for(unsigned int y = 0; y < adjacentVerts[x].size(); y++){
      int q = adjacentVerts[x][y];
      const VEC3F& qpos = vertices[q];
      VEC3F diff = vpos - qpos;
      Real squaredDist = diff.dot(diff) + FORCE_GENERAL_POSITION_EPSILON_SQR;
      double weight = 1.0 / squaredDist;
      sum += weight;
      matrix.insert(x, q) = -weight; 
    }
    matrix.insert(x, x) = sum;
    if(sum <= 0.0001){
      cout << "vertex " << x << "has " << adjacentVerts[x].size() << "surrounding vertices and row has sum " << sum << endl;
      assert(sum > 0);
    }
  }
}

VEC3F* TET_MESH::closestSurfaceNode(const VEC3F& point)
{
  // make sure a surface list was built 
  if (_surfaceVertices.size() == 0) return NULL;

  // tentatively set it to the first in the list
  VEC3F* closest;
  closest = _surfaceVertices[0];
  Real minDistance = (point - *_surfaceVertices[0]).norm();

  // loop through the rest of the vertices
  for (unsigned int x = 1; x < _surfaceVertices.size(); x++)
  {
    // check if this one is closer
    Real distance = (point - *_surfaceVertices[x]).norm();
    if (distance < minDistance)
    {
      minDistance = distance;
      closest = _surfaceVertices[x];
    }
  }
  return closest;
}
void TET_MESH::buildPartitions(vector<int>& tetPartitions)
{
  _tetPartitions = tetPartitions;
  int maxPartition = 0;
  for(unsigned int x = 0; x < tetPartitions.size(); x++){
    if(tetPartitions[x] < 0){
      cout << " some tets have not been assigned a partition!" << endl;
      exit(0);
    }
    maxPartition = (tetPartitions[x] > maxPartition) ? tetPartitions[x] : maxPartition;
  }

  _totalPartitions = maxPartition + 1;

  // tetIDs in each partition
  _partitionedTets.clear();
  _partitionedTets.resize(_totalPartitions);

  for(unsigned int x = 0; x < tetPartitions.size(); x++){
    _partitionedTets[tetPartitions[x]].push_back(x);
  }

  
  _partitionedVertices.clear();
  _partitionedVertices.resize(_totalPartitions);

  vector<vector<pair<int, int> > >originalVertexIDToPartitionIDs(_unconstrainedSize);

  _partitionedOIDToPID.clear();
  _partitionedOIDToPID.resize(_totalPartitions);

  _partitionedSurfaceTetIDs.clear();
  _partitionedSurfaceTetIDs.resize(_totalPartitions);
  _partitionedSurfaceVertexSize.clear();

  _partitionedDofs = 0;
  _partitionDofStartIdx.clear();

  _partitionedSurfaceDofs = 0;
  _partitionSurfaceDofStartIdx.clear();

  for(int p = 0; p < _totalPartitions; p++){
    cout << "partition " << p << endl;
    cout << " number of tets: " << _partitionedTets[p].size() << endl;

    set<int> vertexIDs;
    for(unsigned int x = 0; x < _partitionedTets[p].size(); x++){
      int tetID = _partitionedTets[p][x];
      if(isSurfaceTet(tetID))
        _partitionedSurfaceTetIDs[p].push_back(tetID);

      TET& tet = _tets[tetID];
      for(int y = 0; y < 4; y++){
        int vid = _vertexID[tet.vertices[y]];
        // only care about unconstrained vertices
        if(vid < _unconstrainedSize)
          vertexIDs.insert(vid);
      }
    }

    copy(vertexIDs.begin(), vertexIDs.end(), back_inserter(_partitionedVertices[p]));
    cout << " number of vertices: " << _partitionedVertices[p].size() << endl;

    _partitionDofStartIdx.push_back(_partitionedDofs);
    _partitionedDofs += _partitionedVertices[p].size() * 3;

    int surfaceSize = 0;
    int index = 0;
    while(_partitionedVertices[p][index++] < _surfaceVertexSize)
      surfaceSize++;

    _partitionedSurfaceVertexSize.push_back(surfaceSize);

    _partitionSurfaceDofStartIdx.push_back(_partitionedSurfaceDofs);
    _partitionedSurfaceDofs += surfaceSize * 3;

    for(unsigned int x = 0; x < _partitionedVertices[p].size(); x++){
      int originalVertexID = _partitionedVertices[p][x];

      originalVertexIDToPartitionIDs[originalVertexID].push_back(make_pair(p, x));
      _partitionedOIDToPID[p][originalVertexID] = x;
    }
  }
  _partitionDofStartIdx.push_back(_partitionedDofs);
  _partitionSurfaceDofStartIdx.push_back(_partitionedSurfaceDofs);

  _vertexNumberOfCopies.resize(_unconstrainedSize);
  
  _interfaceSprings.clear();
  _partitionIDs.resize(_unconstrainedSize);

  _surfaceInterfaceSprings.clear();
  _internalInterfaceSprings.clear();

  for(int x = 0; x < _unconstrainedSize; x++){
    vector<pair<int, int> >& partitionIDs = originalVertexIDToPartitionIDs[x];

    _vertexNumberOfCopies[x] = partitionIDs.size();

    int p1 = partitionIDs[0].first;
    int v1 = partitionIDs[0].second;

    _partitionIDs[x] = make_pair(p1, v1);

    for(unsigned int y = 1; y < partitionIDs.size(); y++){
      int p2 = partitionIDs[y].first;
      int v2 = partitionIDs[y].second;
      _interfaceSprings[make_pair(p1, p2)].push_back(make_pair(v1, v2));

      if(x < _surfaceVertexSize){
        _surfaceInterfaceSprings[make_pair(p1, p2)].push_back(make_pair(v1, v2));
      }else{
        _internalInterfaceSprings[make_pair(p1, p2)].push_back(make_pair(v1, v2));
      }
    }
  }

  resetPartitionedAdaptiveMixedSim();
}
void TET_MESH::changeToPartitionOrder(const VECTOR& input, vector<VECTOR>& output)
{
  if(output.size() != _totalPartitions)
    output.resize(_totalPartitions);

  #if USING_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(unsigned int x = 0; x < output.size(); x++){
    output[x].conservativeResize(_partitionedVertices[x].size() * 3);
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      output[x].segment<3>(index * 3) = input.segment<3>(_partitionedVertices[x][y] * 3);
    }
  }
}

void TET_MESH::changeToDefaultPartitionOrder(const VECTOR& input, vector<VECTOR>& output)
{
  if(output.size() != _totalPartitions)
    output.resize(_totalPartitions);

  #if USING_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(unsigned int x = 0; x < output.size(); x++){
    output[x].conservativeResize(_partitionedVertices[x].size() * 3);
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      output[x].segment<3>(y * 3) = input.segment<3>(_partitionedVertices[x][y] * 3);
    }
  }
}

void TET_MESH::changeToPartitionOrder(const VECTOR& input, VECTOR& output)
{
  output.conservativeResize(_partitionedDofs);

  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _totalPartitions; x++){
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      output.segment<3>(_partitionDofStartIdx[x] + index * 3) = input.segment<3>(_partitionedVertices[x][y] * 3);
    }
  }
}
void TET_MESH::restoreDefaultOrder(const VECTOR& input, VECTOR& output)
{
  assert(input.size() == _partitionedDofs);

  output.conservativeResize(_unconstrainedSize * 3);
  output.setZero();

  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _totalPartitions; x++){
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      int vertexID = _partitionedVertices[x][y];
      if(_vertexNumberOfCopies[vertexID] == 1)
        output.segment<3>(vertexID * 3) += input.segment<3>(_partitionDofStartIdx[x] + index * 3);
    }
  }

  for(int x = 0; x < _totalPartitions; x++){
    for(unsigned int y = 0; y < _partitionedVertices[x].size(); y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      int vertexID = _partitionedVertices[x][y];
      if(_vertexNumberOfCopies[vertexID] > 1){
        output.segment<3>(vertexID * 3) += input.segment<3>(_partitionDofStartIdx[x] + index * 3) / _vertexNumberOfCopies[vertexID];
        
      }
    }
  }
}
void TET_MESH::restoreDefaultOrder(const vector<VECTOR>& input, VECTOR& output)
{
  output.conservativeResize(_unconstrainedSize * 3);
  output.setZero();

  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(unsigned int x = 0; x < input.size(); x++){
    for(unsigned int y = 0; y < input[x].size() / 3; y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      int vertexID = _partitionedVertices[x][y];
      if(_vertexNumberOfCopies[vertexID] == 1)
        output.segment<3>(vertexID * 3) += input[x].segment<3>(index * 3);
    }
  }
  for(int x = 0; x < input.size(); x++){
    for(unsigned int y = 0; y < input[x].size() / 3; y++){
      int index = _adaptivePartitionedVertexOrdering[x][y];
      int vertexID = _partitionedVertices[x][y];
      if(_vertexNumberOfCopies[vertexID] > 1){
        output.segment<3>(vertexID * 3) += input[x].segment<3>(index * 3) / _vertexNumberOfCopies[vertexID];
        
      }
    }
  }
}

void TET_MESH::continousExhaustiveCollisionTest(VECTOR& otherDisp, SURFACE* surface)
{
  static Real w = SIMPLE_PARSER::getFloat("continous cd weight", 0.1);

  VECTOR diff = otherDisp.head(_surfaceVertexSize * 3) - _x.head(_surfaceVertexSize * 3);

  for(unsigned int i = 0; i < _surfaceVertices.size(); i++){
    VEC3F position = *_surfaceVertices[i] + w * diff.segment<3>(i * 3);

    if(surface->inside(position)){
      _collisionPairs.push_back(make_pair(_surfaceVertices[i], surface));
    }
  }
  cout << "continous cd: # of collision vertices " << _collisionPairs.size() << endl;
}
void TET_MESH::exhaustiveCollisionTest(SURFACE* surface)
{
  for(unsigned int i = 0; i < _surfaceVertices.size(); i++){
    if(surface->inside(*_surfaceVertices[i]))
      _collisionPairs.push_back(make_pair(_surfaceVertices[i], surface));
  }
  cout << "# of collision vertices " << _collisionPairs.size() << endl;
}
bool TET_MESH::isSurfaceTet(int tetID)
{
  TET& tet = _tets[tetID];
  for(int y = 0; y < 4; y++)
    if(_vertexID[tet.vertices[y]] < _surfaceVertexSize)
      return true;
  return false;
}
void TET_MESH::resetPartitionedAdaptiveMixedSim()
{
  static bool firstCall = true;
  if(firstCall){
    firstCall = false;

    _adaptivePartitionedVertexOrdering.resize(_totalPartitions);
    
    for(int x = 0; x < _totalPartitions; x++)
      _adaptivePartitionedVertexOrdering[x].resize(_partitionedVertices[x].size());

    _partitionedFullsimDofs.resize(_totalPartitions);
    _partitionedReducedsimDofs.resize(_totalPartitions);
    _partitionedFullsimTetIDs.resize(_totalPartitions);
    _partitionFullsimDofStartIdx.resize(_totalPartitions + 1);

    computeVertexAdjacency(_vertexAdjacency);
  }
  
  for(int x = 0; x < _totalPartitions; x++){
    _partitionedFullsimDofs[x] = 0;
    _partitionedReducedsimDofs[x] = _partitionedVertices[x].size() * 3;
    _partitionFullsimDofStartIdx[x] = 0;
    for(int y = 0; y < _partitionedVertices[x].size(); y++)
      _adaptivePartitionedVertexOrdering[x][y] = y;
  }
  _partitionFullsimDofStartIdx[_totalPartitions] = 0;
}

void TET_MESH::initPartitionedHybridSim(vector<int>& contactVertices, VECTOR& isFulsim)
{
  isFulsim.setZero();

  if(contactVertices.size() == 0){
    resetPartitionedAdaptiveMixedSim();
    return;
  }

  vector<bool> visited(_unconstrainedSize, false);
  vector<bool> isActive(_unconstrainedSize, false);

  static Real internalAffectRadius = SIMPLE_PARSER::getFloat("internal active region radius", 0.03);
  static Real surfaceAffectRadius = SIMPLE_PARSER::getFloat("surface active region radius", 0.05);

  queue<int> Q;
  queue<VEC3F> centerQ;

  for(unsigned int x = 0; x < contactVertices.size(); x++){
    Q.push(contactVertices[x]);
    centerQ.push(_restPose[contactVertices[x]]);
  }
  while(!Q.empty()){
    int vertexID = Q.front();
    Q.pop();
    VEC3F center = centerQ.front();
    centerQ.pop();
    if(visited[vertexID])
      continue;

    if((vertexID < _surfaceVertexSize && (_restPose[vertexID] - center).norm() <= surfaceAffectRadius) || (_restPose[vertexID] - center).norm() <= internalAffectRadius){
      isActive[vertexID] = true;
      isFulsim[vertexID] = 1.0;

      vector<int>& neighbors = _vertexAdjacency[vertexID];
      for(unsigned int y = 0; y < neighbors.size(); y++){
        if(neighbors[y] < _unconstrainedSize && !visited[neighbors[y]]){
          Q.push(neighbors[y]);
          centerQ.push(center);
        }
      }
    }
    visited[vertexID] = true;
  }

  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _totalPartitions; x++){
    vector<int>& vertices = _partitionedVertices[x];
    for(unsigned int y = 0; y < vertices.size(); y++){
      _adaptivePartitionedVertexOrdering[x][y] = -1;
    }
    int index = 0;
    set<int> fullsimTetIDs;
    for(unsigned int y = 0; y < vertices.size(); y++){
      if(isActive[vertices[y]])
      {
        _adaptivePartitionedVertexOrdering[x][y] = index++;
        fullsimTetIDs.insert(_tetMembership[vertices[y]].begin(), _tetMembership[vertices[y]].end());
      }
    }

    _partitionedFullsimTetIDs[x].clear();
    for(set<int>::iterator iter = fullsimTetIDs.begin(); iter != fullsimTetIDs.end(); iter++){
      if(_tetPartitions[*iter] == x)
        _partitionedFullsimTetIDs[x].push_back(*iter);
    }

    _partitionedFullsimDofs[x] = index * 3;

    _partitionedReducedsimDofs[x] = (vertices.size() - index) * 3;
    for(unsigned int y = 0; y < vertices.size(); y++){
      if(_adaptivePartitionedVertexOrdering[x][y] == -1)
        _adaptivePartitionedVertexOrdering[x][y] = index++;
    }
    assert(index == vertices.size());
  }
  for(int x = 1; x <= _totalPartitions; x++){
    if(_partitionedReducedsimDofs[x - 1] == 0)
      cout << "partition " << x - 1 << " entirely in full sim" << endl;
    else if(_partitionedFullsimDofs[x - 1] > 0){
      cout << "full sim dofs for partition " << x - 1 << " "<< _partitionedFullsimDofs[x - 1] << endl;
    }
    _partitionFullsimDofStartIdx[x] = _partitionFullsimDofStartIdx[x - 1] + _partitionedFullsimDofs[x - 1];
  }
}

void TET_MESH::initPartitionedHybridSim(vector<int>& contactVertices, VECTOR& keep, VECTOR& isFulsim)
{
  isFulsim.setZero();

  if(contactVertices.size() == 0 && keep.sum() < 1){
    resetPartitionedAdaptiveMixedSim();
    return;
  }

  vector<bool> visited(_unconstrainedSize, false);
  vector<bool> isActive(_unconstrainedSize, false);

  static Real internalAffectRadius = SIMPLE_PARSER::getFloat("internal active region radius", 0.03);
  static Real surfaceAffectRadius = SIMPLE_PARSER::getFloat("surface active region radius", 0.05);

  queue<int> Q;
  queue<VEC3F> centerQ;

  for(unsigned int x = 0; x < contactVertices.size(); x++){
    Q.push(contactVertices[x]);
    centerQ.push(_restPose[contactVertices[x]]);
  }
  while(!Q.empty()){
    int vertexID = Q.front();
    Q.pop();
    VEC3F center = centerQ.front();
    centerQ.pop();
    if(visited[vertexID])
      continue;
    
    if((vertexID < _surfaceVertexSize && (_restPose[vertexID] - center).norm() <= surfaceAffectRadius) || (_restPose[vertexID] - center).norm() <= internalAffectRadius){

      isActive[vertexID] = true;

      vector<int>& neighbors = _vertexAdjacency[vertexID];
      for(unsigned int y = 0; y < neighbors.size(); y++){
        if(neighbors[y] < _unconstrainedSize && !visited[neighbors[y]]){
          Q.push(neighbors[y]);
          centerQ.push(center);
        }
      }
    }
    visited[vertexID] = true;
  }


  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _totalPartitions; x++){
    vector<int>& vertices = _partitionedVertices[x];
    for(unsigned int y = 0; y < vertices.size(); y++){
      _adaptivePartitionedVertexOrdering[x][y] = -1;
    }
    int index = 0;
    set<int> fullsimTetIDs;
    for(unsigned int y = 0; y < vertices.size(); y++){
      if(isActive[vertices[y]] || keep[vertices[y]] > 0.5)
      {
        isFulsim[vertices[y]] = 1.0;
        
        _adaptivePartitionedVertexOrdering[x][y] = index++;
        fullsimTetIDs.insert(_tetMembership[vertices[y]].begin(), _tetMembership[vertices[y]].end());
      }
    }

    _partitionedFullsimTetIDs[x].clear();
    for(set<int>::iterator iter = fullsimTetIDs.begin(); iter != fullsimTetIDs.end(); iter++){
      if(_tetPartitions[*iter] == x)
        _partitionedFullsimTetIDs[x].push_back(*iter);
    }

    _partitionedFullsimDofs[x] = index * 3;

    _partitionedReducedsimDofs[x] = (vertices.size() - index) * 3;
    for(unsigned int y = 0; y < vertices.size(); y++){
      if(_adaptivePartitionedVertexOrdering[x][y] == -1)
        _adaptivePartitionedVertexOrdering[x][y] = index++;
    }
    assert(index == vertices.size());
  }
  for(int x = 1; x <= _totalPartitions; x++){
    cout << "full sim dofs for partition " << x - 1 << " "<< _partitionedFullsimDofs[x - 1] << endl;
    _partitionFullsimDofStartIdx[x] = _partitionFullsimDofStartIdx[x - 1] + _partitionedFullsimDofs[x - 1];
  }
}

void TET_MESH::computeAndWritePartitionedFullsimTetSurfaces(const string& filename)
{
  if(_partitionFullsimDofStartIdx[_totalPartitions] == 0)
    return;

  vector<bool> isFullsimTet(_tets.size(), false);
  for(int x = 0; x < _totalPartitions; x++){
    if(_partitionedFullsimDofs[x] == 0)
      continue;
    vector<int>& tetIDs = _partitionedFullsimTetIDs[x];
    for(unsigned int y = 0; y < tetIDs.size(); y++)
      isFullsimTet[tetIDs[y]] = true;
  }

  // hash the faces according to the sum of their addresses
  map<long long, vector<TRIANGLE> > faceHash;

  // hash the tet and triangle index as well so that everything can be dumped
  // to a file after
  map<long long, vector<int> > tetHash;
  map<long long, vector<int> > triangleHash;
  
  // track the addresses generated, make them long longs to support
  // 64 bit architectures
  vector<long long> sums;

  // insert all the tet faces into the face hash
  for (unsigned int x = 0; x < _tets.size(); x++){
    if(!isFullsimTet[x])
      continue;
    for (int y = 0; y < 4; y++)
    {
      TRIANGLE face = _tets[x].face(y);

      // sum the addresses for use as the hash
      long long sum = 0;
      for (int z = 0; z < 3; z++)
        sum = sum + (long long)(face.vertices[z]);

      // hash it
      faceHash[sum].push_back(_tets[x].face(y));
      tetHash[sum].push_back(x);
      triangleHash[sum].push_back(y);
    }
  }

  vector<TRIANGLE> surfaceFaces;

  // go through all the triangles
  // if more than one hashed in, check for duplicates
  map<long long, vector<TRIANGLE> >::iterator iter;
  vector<int> faceBelongsToTet;
  vector<int> whichFaceInTet;
  int hashSize = faceHash.size();
  int i = 0;
  for (iter = faceHash.begin(); iter != faceHash.end(); iter++, i++)
  {
    vector<TRIANGLE> faces = iter->second;
    if (faces.size() == 1)
    {
      surfaceFaces.push_back(faces[0]);
      faceBelongsToTet.push_back(tetHash[iter->first][0]);
      whichFaceInTet.push_back(triangleHash[iter->first][0]);
      continue;
    }

    // see if this face matches any other ones
    for (unsigned int x = 0; x < faces.size(); x++)
    {
      bool match = false;
      for (unsigned int y = 0; y < faces.size(); y++)
      {
        if (y == x) continue;

        // exploit overloaded TRIANGLE operator
        if (faces[x] == faces[y]) {
          match = true;
          continue;
        }
      }

      // if there are no matches, it is a surface triangle
      if (!match) 
      {
        surfaceFaces.push_back(faces[x]);
        faceBelongsToTet.push_back(tetHash[iter->first][x]);
        whichFaceInTet.push_back(triangleHash[iter->first][x]);
      }
    }
  }

  set<VEC3F*> surfaceVertexSet;
  for(unsigned int x = 0; x < surfaceFaces.size(); x++){
    TRIANGLE& face = surfaceFaces[x];
    for(int y = 0; y < 3; y++){
      surfaceVertexSet.insert(face.vertices[y]);
    }
  }
  vector<VEC3F*> surfaceVertices;
  std::copy(surfaceVertexSet.begin(), surfaceVertexSet.end(), back_inserter(surfaceVertices));
  map<VEC3F*, int> surfaceVertexID;

  for(unsigned int x = 0; x < surfaceVertices.size(); x++){
    surfaceVertexID[surfaceVertices[x]] = x + 1;
  }

  FILE* file = fopen(filename.c_str(), "w");
  if(file == NULL){
    cout << "Cannot open " << filename << " to write!!" << endl;
    return;
  }
  for(int x = 0; x < surfaceVertices.size(); x++){
    VEC3F vertex = *surfaceVertices[x];
    fprintf(file, "v %f %f %f\n", vertex[0], vertex[1], vertex[2]);
  }


  for(int x = 0; x < surfaceFaces.size(); x++){
    TRIANGLE& face = surfaceFaces[x];
    int indices[] = {surfaceVertexID[face.vertices[0]],
                     surfaceVertexID[face.vertices[1]],
                     surfaceVertexID[face.vertices[2]]};
    fprintf(file, "f %d %d %d\n", indices[0], indices[1], indices[2]);
  }
  fclose(file);
}
