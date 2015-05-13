#ifndef SELF_COLLISION_DETECTOR_H
#define SELF_COLLISION_DETECTOR_H

#include <iostream>
#include <SETTINGS.h>
#include <geometry/TET_MESH.h>
#include <geometry/SKELETON.h>
#include <geometry/RIGGER.h>
#include <dtgrid/SPARSE_SDF.h>
#include <geometry/VERTEX_BVH.h>
#include <geometry/TET_BVH.h>
#include <cubature/SCF_CUBATURE_LOADER.h>

using namespace::std;

class SELF_COLLISION_INFO{
public:
  int vertexID;
  int vertexPartition;
  int trianglePartition;
  Real avgArea;
  VEC3I triangleVertexIDs;
  int faceID;

  VEC3F baryCenter;
  bool isCubaturePair;
  MATRIX3 M;
  Real cubatureWeight;
};

class EX_COLLISION_INFO{
public:
  VEC3F penetratingPosition;
  VEC3I triangleVertexIDs;
  VEC3F baryCenter;
  Real avgArea;
  MATRIX3 M;
};

template<class BONE>
class SELF_COLLISION_DETECTOR
{
public:
  SELF_COLLISION_DETECTOR(TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~SELF_COLLISION_DETECTOR();
  void init();

  vector<SELF_COLLISION_INFO>& selfCollisionPoints() { return _selfCollisionPoints; }

  vector<EX_COLLISION_INFO>& externalCollisionPoints() { return _externalCollisionPoints; }

  Real vertexVsTetSCD();

  void drawSelfCollisionPoints();
  void drawExternalCollisionPoints();

  Real collide(vector<VEC3F>& collidingNodes);

private:
  Real findNearestSurfacePoint(int vertexPID, VEC3F penetratingRestPosition, SELF_COLLISION_INFO& info, bool checkProjection = false);

  Real findNearestSurfacePoint(const VEC3F& penetratingRestPosition, EX_COLLISION_INFO& info);

  Real multiDomainVertexVsTetSCD();
  Real singleDomainVertexVsTetSCD();
  Real lowresCubatureSelfCollisionTest();

  void readSCDList();

private:
  TET_MESH* _tetMesh;
  RIGGER<BONE>* _rigger;
  SKELETON<BONE>* _skeleton;

  SPARSE_SDF& _restSDF;
  VERTEX_BVH* _surfaceVertexBVH;
  TET_BVH* _tetBVH;

  vector<VERTEX_BVH*> _partitionedSurfaceVertexBVHs;
  vector<TET_BVH*> _partitionedTetBVHs;

  bool _pseudoCCD;
  bool _useLowresTets;
  bool _useScfCubature;

  SCF_CUBATURE_LOADER* _selfCollisionCubatureLoader;

  map<pair<int, int>, pair<VEC3F, QUATERNION> > _relativeTransforms;

  vector<SELF_COLLISION_INFO> _selfCollisionPoints;
  vector<EX_COLLISION_INFO> _externalCollisionPoints;

  vector<vector<SELF_COLLISION_INFO> > _selfCollisionPointsCopies;

  vector<bool> _selfcheck;
  vector<bool> _ignoreDomain;
  set<pair<int, int> > _ignoreDomainPair;

};

#include "SELF_COLLISION_DETECTOR.inl"

#endif
