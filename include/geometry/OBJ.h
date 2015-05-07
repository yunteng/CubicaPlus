#ifndef _OBJ_H_
#define _OBJ_H_

#include <string>
#include <vector>
#include <map>
#include <SETTINGS.h>
#include <geometry/TRIANGLE.h>
#include <util/MATRIX_UTIL.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
// ===
// OBJ file loader
// ===
// This class can be used to read a .obj 3d model.  
// 
// The model is stored as a sequence of faces with vertices and normals.
// 
// Each face contains an array of idices into the vertices and normals of the
// OBJ.  These arrays are always to be of the same size.
// 
// The number of points in the ith faces is
//  faces[i].vertices.size ()
// 
// For the jth point (index starting at 0) in face i (index starting at 0), its position, normal
// and texture coordinate are respectively
//  vertices[faces[i].vertex[j]]
//  normals[faces[i].normals[j]]
//  texcoords[faces[i].texcoords[j]]
// 
// However, if the normal or texcoord of a vertex is not specified in the OBJ file, then the
// value of faces[i].normals[j] or faces[i].texcoords[j] will be -1 (and the above will not hold)
//////////////////////////////////////////////////////////////////////////////////
class OBJ
{
public:
  // A single face
  // The vertices/texcoords/normals arrays are all of the same size and contain 
  // indices into the OBJ's vertices/normals/texcoords arrays.  If the normal or 
  // texcoord of a point on a face is not specified by the OBJ file then its
  // index will be -1
  struct Face
  {
    std::vector<int> vertices;
    std::vector<int> texcoords;
    std::vector<int> normals;

    void appendVertFrom(int vSrc, Face& src)
    {
      vertices.push_back(src.vertices[ vSrc ]);
      texcoords.push_back(src.texcoords[ vSrc ]);
      normals.push_back(src.normals[ vSrc ]);
    } 
  };

  // return true if succeeded to read obj from fileName
  bool Load(const std::string& fileName, Real scale = 1.0);

  // return true if succeeded to write obj to fileName
  bool Save(const std::string& fileName);
  
  bool SaveToPBRT (const string& fileName);

  // return true if succeeded to write obj to fileName
  bool SaveFiltered(vector<int>& filteredFaces, const std::string& fileName);

  // return the max/min coordinates of a bounding box around the object
  void BoundingBox (VEC3F& minPosition, VEC3F& maxPosition);
  void BoundingBox(Real* boundingBox) {
    VEC3F minVertex;
    VEC3F maxVertex;
    BoundingBox(minVertex, maxVertex);
    boundingBox[0] = minVertex[0]; boundingBox[1] = maxVertex[0];
    boundingBox[2] = minVertex[1]; boundingBox[3] = maxVertex[1];
    boundingBox[4] = minVertex[2]; boundingBox[5] = maxVertex[2];
  };

  // compute some reasonable vertex normals
  void ComputeVertexNormals();
  void SmoothVertexNormals();
  void FlipVertexNormals();

  // return true if succeeded to read PLY from fileName
  // bool LoadPly(const std::string& fileName);

  // return true if read in a PBRT file
  bool LoadPBRT(string& filename);

  // return true if succeeded to read PLY from fileName
  // bool SavePly(const std::string& fileName);

  // Returns the i-th vertex of the given face (assuming it's from this object)
  VEC3F GetVertex(const Face& f, unsigned int i) const;

  // Get various face statistics
  VEC3F GetCentroid(const Face& f) const;
  Face& ClosestFace(const VEC3F& p);
  float DistanceToClosestFace(const VEC3F& p);
  int ClosestVertex(const VEC3F& p);

  // vertex positions
  std::vector<VEC3F> vertices;

  // vertex normals
  std::vector<VEC3F> normals;

  // vertex texture coordinates
  std::vector<VEC2F> texcoords;

  // mesh faces
  std::vector<Face> faces;

  // get a specific vertex on a face
  VEC3F getFaceVert(int iface, int ifaceVert)
  {
    return vertices[ faces[iface].vertices[ifaceVert] ];
  }

  // draw everything to GL
  void drawVert(int i);
  void draw();
  void drawFiltered();
  void drawFilteredSubset();

  // fit the geometry inside a 1x1x1 box
  void normalize(int res = 0);

  // apply the same normalization to this vertex -- assumes normalize(res) has been called!!!
  void normalize(VEC3F& vertex);

  // tell the OBJ what the res of the BCC grid is
  void setBCCRes(int bccRes) { _bccRes = bccRes;};

  // inside / outside tests
  bool insideExhaustive(float* point);

  float distance(float* point);
  float bruteForceDistance(float* point);

  bool intersect_triangle(float orig[3], float dir[3], int faceIndex);

  // translate the whole mesh
  void translate(VEC3F translation);

  void createAccelGrid();
  void createOccupancyGrid(int fatten = 2);
  void createDistanceGrid(int res);

  const string& filename() { return _filename; };
  Real pointFaceDistanceSq(Face& f, const VEC3F& point);

  std::vector<int>& vertexGroups() { return _vertexGroupStarts; };
  std::vector<int>& faceGroups() { return _faceGroupStarts; };
  std::vector<int>& texcoordGroups() { return _texcoordGroupStarts; };
  std::vector<string>& groupNames() { return _groupNames; };
  std::vector<int>& materialStarts() { return _materialStarts; };

  OBJ() : _filename("") { _scaling = -1.0; _bccRes = -1; };
  OBJ(const OBJ& toCopy);

  VEC3F center() {
    VEC3F c;
    for (unsigned int x = 0; x < vertices.size(); x++)
      c += vertices[x];
    c *= 1.0 / vertices.size();
    return c;
  };

  void transformMesh(MATRIX3 transform);
  void translateMesh(VEC3F translation);
  
  void scale(Real s);

  void drawAccelGrid();
  void drawOccupancyGrid();

  // threshold at which to cull stretched triangles
  float& filteringThreshold() { return _filteringThreshold; };

  std::vector<float>& triangleAreas()  { return _triangleAreas; };
  std::vector<float>& maxEdgeLengths() { return _maxEdgeLengths; };

  // load in a TET_MESH surface mesh, just to fake out the AOV renderer
  void LoadTetSurfaceMesh(vector<TRIANGLE*> faces);

  // refresh a TET_MESH surface mesh, just to fake out the AOV renderer
  void RefreshTetSurfaceMesh(vector<TRIANGLE*> faces);

private:
  // PLY types, only used in LoadPly
  // struct PlyVertex {
  //   float x,y,z;             // the usual 3-space position of a vertex
  // };

  // struct PlyFace {
  //   unsigned char intensity; // this user attaches intensity to faces
  //   unsigned char nverts;    // number of vertex indices in list
  //   int *verts;              // vertex index list
  // };

  // resolution of the isostuffer -- used to set the resolution of 
  // the distance and accel grids
  int _bccRes;

  // a grid of vectors -- given a grid index, return the vector
  // of faces whose AABB intersects that grid cell. Used to speed up
  // inside/outside test
  map<int, vector<int> > _accelGrid;

  // grid if distance values -- given a grid index, return the distance
  // to the closest face
  map<int, Real> _distanceGrid;

  // grid of occupany values -- given a grid index, return a bool
  // that tells if a triangle overlapped this grid cell
  map<int, bool > _occupancyGrid;

  // use the accelGrid to lookup the closest face
  Real fastDistToClosestFace(int vertexID);

  // dimensions of the accelGrid
  int _xRes, _yRes, _zRes;

  // dimensions of the distance grid
  int _xResDist, _yResDist, _zResDist;

  string _filename;

  std::vector<int> _vertexGroupStarts;
  std::vector<int> _faceGroupStarts;
  std::vector<int> _texcoordGroupStarts;
  std::vector<string> _groupNames;
  std::vector<int> _materialStarts;

  // store normalization params
  VEC3F _centerOfMass;
  Real  _scaling;

  // cache original face areas and maximum edge lengths --
  // if the embedding of a mesh is bad, the triangles will stretch too much, and
  // they should be filtered out
  std::vector<float> _triangleAreas;
  std::vector<float> _maxEdgeLengths;
  float _filteringThreshold;
};

#endif
