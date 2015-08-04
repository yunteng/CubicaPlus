/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <dtgrid/SPARSE_SDF.h>

int SPARSE_SDF::trianglesAt(const VEC3F& position, TRIANGLE** &triangles)
{
  MyLocator loc;
  MyVec3 pos(position[0], position[1], position[2]);
  if(_grid->getLocator(pos, &loc)){
      triangles = _hashSurface[loc.iv3D];
  }
  return _gridTriangleNum[loc.iv3D];
}
SPARSE_SDF::~SPARSE_SDF()
{
  for(unsigned int x = 0; x < _hashSurface.size(); x++){
    if(_hashSurface[x] != NULL){
      delete[] _hashSurface[x];
    }
  }
  delete _grid;
}
void SPARSE_SDF::hashTriangle(TRIANGLE& triangle, map<int, vector<TRIANGLE*> >& tempHash)
{
  VEC3F mins, maxs;
  triangle.boundingBox(mins, maxs);
  MyVec3 dtMin(mins[0], mins[1], mins[2]);
  MyVec3 dtMax(maxs[0], maxs[1], maxs[2]);

  Index iMins[3];
  Index iMaxs[3];
  _grid->boundingBoxIndices(dtMin, dtMax, iMins, iMaxs);

  for(Index z = iMins[2]; z <= iMaxs[2]; z++)
    for (Index y = iMins[1]; y <= iMaxs[1]; y++)
      for (Index x = iMins[0]; x <= iMaxs[0]; x++){
        MyLocator loc;
        bool gridExist = _grid->getLocator(x, y, z, &loc);
        if(gridExist)
            tempHash[loc.iv3D].push_back(&triangle);
      }
}

void SPARSE_SDF::initHashSurface(vector<TRIANGLE>& surfaces)
{
  _hashSurface.resize(_grid->getNumVa3D(), NULL);
  _gridTriangleNum.resize(_grid->getNumVa3D(), 0);

  map<int, vector<TRIANGLE*> > tempHash;
  for(int i = 0; i < surfaces.size(); i++){
    hashTriangle(surfaces[i], tempHash);
  }
  for(map<int, vector<TRIANGLE*> >::iterator iter = tempHash.begin(); iter != tempHash.end(); iter++){
    int index = iter->first;
    vector<TRIANGLE*>& triangles = iter->second;

    _hashSurface[index] = new TRIANGLE*[triangles.size()];
    _gridTriangleNum[index] = triangles.size();
    for(unsigned int x = 0; x < triangles.size(); x++){
      _hashSurface[index][x] = triangles[x];
    }
  }
}

bool SPARSE_SDF::nearestSurfacePoint(const VEC3F& position, VEC3F& surface)
{
    MyVec3 vecNormal;
    Real sdf;

    MyVec3 positionVec(position[0], position[1], position[2]);
    MyVec3 tmpSurface = positionVec;

    bool inNarrowBand = _grid->lookupDistanceAndNormalFast(positionVec, sdf, vecNormal);
    if(sdf > 0){
        return false;
    }
    if(!inNarrowBand){
        cout << "lookup point inside the mesh but not captured by the narrow band!!!" << endl;
        // Core::throwDefaultException(std::string("lookup point inside the mesh but not captured by the narrow band!!!"), __FILE__, __LINE__);
        return false;
    }

    Real initSdf = sdf;
    int iter = 0;

    while(abs(sdf) > 1e-9 && iter < 100){
        tmpSurface +=  (-sdf) * vecNormal;
        sdf = (*_grid)(tmpSurface);
        iter++;
    }
  
    Real dist = (tmpSurface - positionVec).length();
    if(dist > abs(initSdf * 1.5)){
        return false;
    }
    surface[0] = tmpSurface[0];
    surface[1] = tmpSurface[1];
    surface[2] = tmpSurface[2];
    return true;
}

bool SPARSE_SDF::absoluteNearestSurfacePoint(const VEC3F& position, VEC3F& surface)
{
    MyVec3 vecNormal;
    Real sdf;

    MyVec3 positionVec(position[0], position[1], position[2]);
    MyVec3 tmpSurface = positionVec;

    bool inNarrowBand = _grid->lookupDistanceAndNormalFast(positionVec, sdf, vecNormal);

    if(!inNarrowBand){
        cout << "lookup point inside the mesh but not captured by the narrow band!!!" << endl;
        // Core::throwDefaultException(std::string("lookup point inside the mesh but not captured by the narrow band!!!"), __FILE__, __LINE__);
        return false;
    }

    Real initSdf = sdf;
    int iter = 0;

    while(abs(sdf) > 1e-9 && iter < 100){
        tmpSurface +=  (-sdf) * vecNormal;
        sdf = (*_grid)(tmpSurface);
        iter++;
    }
  
    Real dist = (tmpSurface - positionVec).length();
    if(dist > abs(initSdf * 1.5)){
        return false;
    }
    surface[0] = tmpSurface[0];
    surface[1] = tmpSurface[1];
    surface[2] = tmpSurface[2];
    return true;
}
