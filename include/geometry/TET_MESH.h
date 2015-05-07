#ifndef TET_MESH_H
#define TET_MESH_H

#include <SETTINGS.h>
#include <util/SIMPLE_PARSER.h>
#include <geometry/TRIANGLE.h>
#include <geometry/SURFACE.h>
#include <geometry/TET.h>
#include <dtgrid/SPARSE_SDF.h>
#include <linearalgebra/COO_MATRIX.h>
#include <linearalgebra/BLOCK_COO_MATRIX.h>
#include <linearalgebra/BLOCK_SPARSE_MATRIX.h>
#include <fstream>

#if USING_OPENMP
#include <omp.h>
#endif

using namespace::std;

class TET_MESH
{
public:
  TET_MESH();

  /*
  initialize functions
  */
  void init();
  void computeSurfaceFaces();
  void computeSurfaceVertexNormal();
  bool readSurfaceFaceCache();

  void computeMasses();
  inline COO_MATRIX& massMatrix() { return _masses; };
  inline VECTOR& massVector() { return _massVec; };
  inline Real mass(int vertexID) { assert(vertexID < _massVec.size()); return _massVec[vertexID]; };
  
  void drawSurfaceFaces();
  void drawPartition(int partition);
  void drawConstrainedNodes();
  void drawLowresEmbedding();
  void drawDisplacement();
  void drawSurfaceVertexOneRing(int surfaceVertexID);
  void drawCollisionPairs();
  void drawInterfaceVertices();
  void drawFullsimVertices();
  void click(const VEC3F& point)
  {
    _clickedNode = closestSurfaceNode(point);
  }
  void unclick() { _clickedNode = NULL; };
  

  /*
  accesses
  */
  inline const string& filename() { return _filename; }; 
  inline vector<VEC3F>& restPose() { return _restPose; };
  inline vector<VEC3F>& vertices() { return _vertices; };
  inline vector<TET>& tets() { return _tets; };
  inline vector<TRIANGLE>& surfaceFaces() { return  _surfaceFaces; };
  inline vector<VEC3F*>& surfaceVertices() { return _surfaceVertices; };
  inline Real restSurfaceFaceArea(TRIANGLE* triangle) { 
    map<TRIANGLE*, Real>::iterator ptr = _restSurfaceFaceAreas.find(triangle);
    if(ptr == _restSurfaceFaceAreas.end()){
      cout << "Look up triangle is not a surface face!!" << endl;
      return -1;
    }
    return ptr->second;
  };
  inline Real restSurfaceVertexArea(VEC3F* vertex) { map<VEC3F*, Real>::iterator ptr = _restSurfaceVertexAreas.find(vertex);
    if(ptr == _restSurfaceVertexAreas.end()){
      cout << "Look up vertex is not on the surface!!";
      return -1;
    }
    return ptr->second;
  };

  int totalTets() const          { return (int)(_tets.size()); };

  inline int surfaceVertexSize()   { return _surfaceVertexSize; };
  inline int unconstrainedNodes()  { return _unconstrainedSize; };
  inline bool constrained() { return _constrainedSize > 0; };
  inline bool isConstrained(VEC3F* vertex) { return _vertexID[vertex] >= _unconstrainedSize; };
  inline bool isConstrained(int vid) { return vid >= _unconstrainedSize; };
  inline bool isSurfaceVertex(int vid) { return vid < _surfaceVertexSize; };
  inline bool isSurfaceVertex(VEC3F* vertex) { return _vertexID[vertex] < _surfaceVertexSize; };
  
  inline int totalCores() { return _totalCores;}

  inline int dofs() { return _unconstrainedSize * 3; };
  inline VECTOR& x() { return _x; };
  inline vector<MATRIX3>& F() { return _F; };
  inline VECTOR& internalForce() { return _R; };
  inline VECTOR& surfaceInternalForce() { return _surfaceR; }
  inline COO_MATRIX& stiffnessMatrix() { return _stiffness; };
  inline BLOCK_SPARSE_MATRIX& surfaceStiffnessMatrix() { return _surfaceStiffness;}

  int& totalMaterials()          { return _totalMaterials; };
  MATERIAL*** materialCopies() { return _materialCopies; };

  int vertexID(VEC3F* vertex)      { return (_vertexID.find(vertex) != _vertexID.end()) ? _vertexID[vertex] : -1; };
  int surfaceFaceID(TRIANGLE* face) { return _surfaceFaceID[face]; };

  VEC3F* vertex(int id)     { return &_vertices[id]; };
  VEC3F* restVertex(int id) { return &_restPose[id]; };

  SPARSE_SDF& restSDF()     { return _restSDF; };
  /*
  IO function
  */
  bool readTetGen(string prefix);
  bool readMeshFile(const string& filename);
  void readMeshFile(FILE* file);
  void writeMeshFile(const string& filename);
  void writeMeshFile(FILE* file);
  void writeObj(const string& filename);
  void writeNewConstraints(vector<VEC3F*>& newConstraints);
  void readMaterialsIndices();
  void writeMaterialIndices();
  void loadMaterials();

  void readDisplacementFromRest(const string& filename);
  void writeDisplacementFromRest(const string& filename);
  void writeSurfaceVertexPositions(const string& filename);
  void getSurfaceVertexDisp(VECTOR& sv);

  /*
  use a low res embedded mesh to accelerate scd
  */
  bool readLowresEmbeddedMesh(const string& filename);
  void writeLowresEmbeddedMesh(const string& filename);
  bool readLowresEmbeddingMap(const string& filename);
  void writeLowresEmbeddingMap(const string& filename);
  void mapToLowresEmbedding();

  inline vector<VEC3F>& lowresVertices() { return _lowresVertices; };
  inline vector<VEC3F>& lowresRestPose() { return _lowresRestPose; };
  inline vector<TET>&   lowresTets()     { return _lowresTets;     };
  inline vector<int>&   lowToHighID()    { return _lowToHighID;    };
  int lowresVertexID(VEC3F* vertex)      { return (_lowresVertexID.find(vertex) != _lowresVertexID.end()) ? _lowresVertexID[vertex] : -1; };
  void updateLowresEmbedding();

  void exhaustiveCollisionTest(SURFACE* surface);
  void continousExhaustiveCollisionTest(VECTOR& otherDisp, SURFACE* surface);
  vector<pair<VEC3F*, SURFACE*> >& collisionPairs() { return _collisionPairs; };

  /*
  geometry related
  */
  VEC3F* closestSurfaceNode(const VEC3F& point);
  VEC3F* clickedNode() { return _clickedNode; };
  void smoothSurface();

  void computeSurfaceVertexOneRings();

  void computeVertexAdjacency(vector<vector<int> >& adjacentVerts);
  void computeLaplacianMatrix(const vector<vector<int> >& adjacentVerts, SpMat& matrix);

  /*
  simualtion related
  */
  void generateF();

  void recoverX();
  void updateFullMesh();
  void updateFullMesh(const VECTOR& pos);

  void reset();

  /*
  partitioning
  */
  void buildPartitions(vector<int>& tetPartitions);
  inline int totalPartitions() { return _totalPartitions; };

  inline vector<vector<int> >& partitionedTets() { return _partitionedTets; };

  inline vector<int>& partitionedTets(int partition){ return _partitionedTets[partition]; };

  inline map<int, int>& partitionedVertexIDs(int partition) { return _partitionedOIDToPID[partition]; };

  inline int partitionedVertexID(int partition, int vertexID) { 
    if(_partitionedOIDToPID[partition].find(vertexID) == _partitionedOIDToPID[partition].end())
      return -1;
    return _partitionedOIDToPID[partition][vertexID]; 
  };
  inline pair<int, int>& partitionedVertexID(int vertexID) { return _partitionIDs[vertexID]; };

  inline map<pair<int, int>, vector<pair<int, int> > >& interfaceSprings() { return _interfaceSprings; };

  inline map<pair<int, int>, vector<pair<int, int> > >& internalInterfaceSprings() { return _internalInterfaceSprings; };

  inline map<pair<int, int>, vector<pair<int, int> > >& surfaceInterfaceSprings() { return _surfaceInterfaceSprings; };

  inline vector<vector<int> >& partitionedVertices() { return _partitionedVertices; };

  inline vector<int>& partitionedVertices(int x) { return _partitionedVertices[x]; };
  inline vector<int>& partitionedSurfaceTetIDs(int x) { return _partitionedSurfaceTetIDs[x]; };
  inline int partitionedSurfaceVertexSize(int x) { return _partitionedSurfaceVertexSize[x]; };

  inline int partitionedDofs() { return _partitionedDofs; };
  inline int partitionDofStartIdx(int x) { return _partitionDofStartIdx[x]; };

  inline int partitionedSurfaceDofs() { return _partitionedSurfaceDofs; };
  inline int partitionSurfaceDofStartIdx(int x) { return _partitionSurfaceDofStartIdx[x]; };

  void changeToPartitionOrder(const VECTOR& input, vector<VECTOR>& output);
  void changeToDefaultPartitionOrder(const VECTOR& input, vector<VECTOR>& output);
  void changeToPartitionOrder(const VECTOR& input, VECTOR& output);
  void restoreNatualOrder(const vector<VECTOR>& input, VECTOR& output);
  void restoreNatualOrder(const VECTOR& input, VECTOR& output);

  /*
  surface condensation
  */
  vector<int>& surfaceTetIDs() {return _surfaceTetIDs; }
  bool isSurfaceTet(int tetID);

  inline VECTOR& partitionedInternalForce() { return _partitionedR; };
  inline VECTOR& partitionedSurfaceInternalForce() { return _partitionedSurfaceR; }
  inline BLOCK_COO_MATRIX& partitionedStiffnessMatrix() { return _partitionedStiffness; };
  inline BLOCK_SPARSE_MATRIX& partitionedFullStiffnessDiag() { return _partitionedFullStiffnessDiag; };
  inline BLOCK_SPARSE_MATRIX& partitionedFullStiffnessOffDiag() { return _partitionedFullStiffnessOffDiag; };
  inline BLOCK_SPARSE_MATRIX& partitionedReducedStiffnessDiag() { return _partitionedReducedStiffnessDiag; };

  inline bool isNeighbors(int x, int y) {
    return _interfaceSprings.find(make_pair(x, y)) != _interfaceSprings.end() || _interfaceSprings.find(make_pair(y, x)) != _interfaceSprings.end();
  };

  void resetPartitionedAdaptiveMixedSim();
  void initPartitionedAdaptiveMixedSim(vector<int>& contactVertices, VECTOR& isFulsim);
  void initPartitionedAdaptiveMixedSim(vector<int>& contactVertices, VECTOR& keepFullsim, VECTOR& isFulsim);

  inline vector<int>& adaptivePartitionedVertexOrdering(int x) { return _adaptivePartitionedVertexOrdering[x]; };
  inline vector<int>& partitionedFullsimDofs() { return _partitionedFullsimDofs; };
  inline int partitionedFullsimDofs(int x) { return _partitionedFullsimDofs[x]; };
  inline vector<int>& partitionedReducedsimDofs() { return _partitionedReducedsimDofs; };
  inline int partitionedReducedsimDofs(int x) { return _partitionedReducedsimDofs[x]; };
  inline vector<int>& partitionedFullsimTetIDs(int x) { return _partitionedFullsimTetIDs[x]; };
  inline int partitionFullsimDofStartIdx(int x) { return _partitionFullsimDofStartIdx[x]; };

  inline bool isPartitionedFullsimVertex(int p, int vID)
  {
    return _partitionedFullsimDofs[p] > 0 && _adaptivePartitionedVertexOrdering[p][vID] * 3 < _partitionedFullsimDofs[p];
  }
  inline bool isPartitionedFullsimVertex(const pair<int, int>& partitionID)
  {
    return _partitionedFullsimDofs[partitionID.first] > 0 && _adaptivePartitionedVertexOrdering[partitionID.first][partitionID.second] * 3 < _partitionedFullsimDofs[partitionID.first];
  }
  inline int adaptivePartitionedVertexID(const pair<int, int>& partitionID){
    return _adaptivePartitionedVertexOrdering[partitionID.first][partitionID.second];
  }
  inline bool isCloned(int vertexID) { return _vertexNumberOfCopies[vertexID] > 1; }

  void interpolateToLowresMesh(const VECTOR& input, VECTOR& output);
  void interpolateFromLowresMesh(const VECTOR& input, VECTOR& output);

  void computeAndWritePartitionedFullsimTetSurfaces(const string& filename);

protected:
  string _filename;
  vector<VEC3F> _vertices;
  vector<VEC3F> _restPose;
  VECTOR _surfaceNormals;

  SPARSE_SDF _restSDF;

  COO_MATRIX _masses;
  VECTOR _massVec;
  Real _totalMass;
  Real _totalVolume;

  // _vertices[_unconstrainedSize] to 
  // _vertices[_unconstrainedSize + _constrainedSize - 1], 
  // are constrained
  int _constrainedSize;

  // _vertices[0] to _vertices[_unconstrainedSize - 1]
  // are unconstrained
  int _unconstrainedSize;

  // _vertices[0] to _vertices[_surfaceVertexSize - 1]
  // are surfaceVertices, assume no surface vertex is constrained
  int _surfaceVertexSize;

  // maps vertex address to its index in _vertices
  map<VEC3F*, int> _vertexID;

  map<int, vector<int> > _tetMembership;

  vector<TET> _tets;
  vector<TRIANGLE> _surfaceFaces;
  map<TRIANGLE*, int> _surfaceFaceID;
  vector<VEC3F*> _surfaceVertices;
  vector<bool>   _isSurfaceVertex;
  map<TRIANGLE*, Real> _restSurfaceFaceAreas;
  map<VEC3F*, Real> _restSurfaceVertexAreas;

  vector<vector<VEC3F*> > _surfaceVertexOneRings;

  // lowres embedding, used for self collision detection
  int _lowSurfaceVertexSize;
  vector<VEC3F> _lowresVertices;
  vector<VEC3F> _lowresRestPose;
  vector<TET>   _lowresTets;
  vector<TRIANGLE> _lowresSurfaceFaces;
  vector<int>   _lowToHighID;
  vector<pair<int, QUATERNION> > _highToLowInternalBaryCentricCoordinates;
  vector<pair<int, VEC3F> > _highToLowSurfaceBaryCentricCoordinates;
  
  vector<pair<int, QUATERNION> > _lowToHighInternalBaryCentricCoordinates;
  vector<pair<int, VEC3F> > _lowToHighSurfaceBaryCentricCoordinates;

  map<VEC3F*, int> _lowresVertexID;

  int _totalCores;
  MATERIAL*** _materialCopies;
  int _totalMaterials;

  VEC3F* _clickedNode;

  COO_MATRIX _stiffness;

  BLOCK_SPARSE_MATRIX _surfaceStiffness;

  VECTOR _x;
  vector<MATRIX3> _F;
  VECTOR _R;
  VECTOR _surfaceR;

  vector<pair<VEC3F*, SURFACE*> > _collisionPairs;

  map<VEC3F*, int> _surfaceVertexID;

  // partitioned variables
  int _totalPartitions;
  vector<vector<int> > _partitionedTets;
  vector<vector<int> > _partitionedVertices;
  vector<map<int, int> > _partitionedOIDToPID;
  // for each vertex, just keep one partitionID
  vector<pair<int, int> > _partitionIDs;
  vector<int> _vertexNumberOfCopies;
  map<pair<int, int>, vector<pair<int, int> > > _interfaceSprings;

  map<pair<int, int>, vector<pair<int, int> > > _internalInterfaceSprings;

  map<pair<int, int>, vector<pair<int, int> > > _surfaceInterfaceSprings;


  int _partitionedDofs;
  vector<int> _partitionDofStartIdx;
  vector<int> _tetPartitions;

  /* 
  condensation
  */
  vector<int> _surfaceTetIDs;
  vector<vector<int> > _partitionedSurfaceTetIDs;
  vector<int> _partitionedSurfaceVertexSize;

  int _partitionedSurfaceDofs;
  vector<int> _partitionSurfaceDofStartIdx;

  BLOCK_COO_MATRIX _partitionedStiffness;
  VECTOR     _partitionedR;
  BLOCK_SPARSE_MATRIX _partitionedFullStiffnessDiag;
  BLOCK_SPARSE_MATRIX _partitionedFullStiffnessOffDiag;
  BLOCK_SPARSE_MATRIX _partitionedReducedStiffnessDiag;
  VECTOR     _partitionedSurfaceR;

  // used by partitioned hybrid simulation
  // assign each partititioned vertex a new ID
  // the purpose is to put fullsim vertices in
  // front of reduced sim vertices
  vector<vector<int> > _adaptivePartitionedVertexOrdering;
  vector<int> _partitionedFullsimDofs;
  vector<int> _partitionedReducedsimDofs;
  vector<vector<int> > _partitionedFullsimTetIDs;
  vector<int> _partitionFullsimDofStartIdx;
  Real _maxContactNodeGeodist;

  vector<vector<int> > _vertexAdjacency;
};
#endif
