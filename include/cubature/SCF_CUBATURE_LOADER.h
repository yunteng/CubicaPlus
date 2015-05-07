#ifndef SCF_CUBATURE_LOADER_H
#define SCF_CUBATURE_LOADER_H
#include <algorithm>
#include <queue>
#include <SETTINGS.h>
#include <util/MATRIX_UTIL.h>
#include <util/IO.h>
#include <geometry/TET_MESH.h>

class SCF_CUBATURE
{
public:
  int vertexPartition;
  int trianglePartition;
  int vertexID;
  int surfaceID;
  Real weight;
  void read(FILE* file);
};
// struct orderByValueGreaterThan{
//   bool operator()(pair<int, pair<Real, Real> > const& a, pair<int, pair<Real, Real> > const& b) const{
//     return (a.second.first + a.second.second) > (b.second.first + b.second.second);
//   }
// };

struct orderByValueGreaterThan{
  bool operator()(pair<int, Real > const& a, pair<int, Real > const& b) const{
    return a.second > b.second;
  }
};

class PAIRWISE_SCF_CUBATURES
{
public:
  PAIRWISE_SCF_CUBATURES();
  ~PAIRWISE_SCF_CUBATURES();

  bool blendCubatures(VECTOR& currentCoords, vector<SCF_CUBATURE>& output, bool isCheb = false);

  // vector<vector<SCF_CUBATURE> >& cubatures() { return _cubatures; };

  inline bool& isNeighbor()              { return _isNeighbor; };
  inline Real& safeDistance()            { return _safeDistance; };
  inline int& leftPartition()            { return _leftPartition; };
  inline int& rightPartition()           { return _rightPartition; };

  inline bool cacheValid(const VEC3F& translation, const QUATERNION& rotation) { 
    if((translation - _cachedTranslation).squaredNorm() > 1e-6){
      _cachedTranslation = translation;
      _cachedRotation = rotation;
      return false;
    }

    QUATERNION diff = rotation - _cachedRotation;
    if(diff.x() * diff.x() + diff.y() * diff.y() + diff.z() * diff.z() + diff.w() * diff.w() > 1e-6){
      _cachedRotation = rotation;
      return false;
    }
    return true;
  };
  inline vector<pair<int, Real> >& cachedLeftVertices() { return _cachedLeftVertices; };
  inline vector<pair<int, Real> >& cachedRightVertices() { return _cachedRightVertices; };
  inline bool& cachedFoundMatch() { return _cachedFoundMatch; };
  
  void read(string filename);


  
private:
  MATRIX _pcaBasis;
  vector<VECTOR> _sampleCoordinates;

  vector<vector<SCF_CUBATURE> > _cubatures;

  bool _isNeighbor;
  Real _safeDistance;
  
  int _leftPartition;
  int _rightPartition;


  VEC3F _cachedTranslation;
  QUATERNION _cachedRotation;
  bool _cachedFoundMatch;
  vector<pair<int, Real> > _cachedLeftVertices;
  vector<pair<int, Real> > _cachedRightVertices;
  // vector<Real> _sampleNorms;

  // vector<VEC3F> _translations;
  // vector<QUATERNION> _rotations;
  
};

class SCF_CUBATURE_LOADER
{
public:
  SCF_CUBATURE_LOADER(TET_MESH* tetMesh);
  ~SCF_CUBATURE_LOADER();
  bool loadAllCubatures(const string& dirName);
  
  bool getCubatureSurfaceVertices(pair<int, int> partitionPair, const VEC3F& relativeTranslation, const QUATERNION& relativeRotation, vector<pair<int, Real> >& leftVertices, vector<pair<int, Real> >& rightVertices);

private:
  TET_MESH* _tetMesh;
  bool _isCheb;
  map<pair<int, int>, PAIRWISE_SCF_CUBATURES> _pairwiseCubatures;

  VECTOR getTransformedColumn(int rightPartition, const VEC3F& relativeTranslation, const QUATERNION& relativeRotation);

};
#endif
