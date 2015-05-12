#include <geometry/TET_MESH.h>

TET_MESH::TET_MESH():
  _clickedNode(NULL),
  _totalPartitions(1)
{
#if USING_OPENMP
  _totalCores = omp_get_max_threads();
#else
  _totalCores = 1;
#endif
  cout << "Using " << _totalCores << " cores" << endl;
  string path = SIMPLE_PARSER::getString("output path", "./");
  string meshName = SIMPLE_PARSER::getString("tet mesh name", "");
  string prefix = path + meshName;
  _filename = prefix + ".tetmesh";

  string constrainedMeshFile = _filename + ".constrained";

  if(!readMeshFile(constrainedMeshFile)){
    if(!readMeshFile(_filename)){
      cout << "try reading from tetgen" << endl;
      if(!readTetGen(prefix)){
        exit(0);
      }else{
        writeMeshFile(_filename);
        // reload from the output as the vertex ordering is changed
        readMeshFile(_filename);
      }
    }
  }else
    _filename = constrainedMeshFile;

  init();
}
void TET_MESH::init()
{
  bool simulate = SIMPLE_PARSER::getBool("simulate full", true);
  bool dynamic = SIMPLE_PARSER::getBool("dynamic", false);
  bool verbose = SIMPLE_PARSER::getBool("verbose", false);

  _restPose.resize(_vertices.size());
  for(unsigned int x = 0; x < _vertices.size(); x++)
    _restPose[x] = _vertices[x];

  // record the vertex positions of each vertex
  for(unsigned int x = 0; x < _vertices.size(); x++)
    _vertexID[&_vertices[x]] = x; 

  _tetMembership.clear();
  for(unsigned int x = 0; x < _tets.size(); x++){
    for(int y = 0; y < 4; y++){
      int vid = _vertexID[_tets[x].vertices[y]];
      _tets[x].setRest(y, &_restPose[vid]);
      _tetMembership[vid].push_back(x);
    }
  }
  computeSurfaceFaces();
  _surfaceFaceID.clear();
  for(unsigned int x = 0; x < _surfaceFaces.size(); x++){
    _surfaceFaceID[&_surfaceFaces[x]] = x;
  }

  _surfaceVertexID.clear();
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++){
    _surfaceVertexID[_surfaceVertices[x]] = x;
    assert(_vertexID[_surfaceVertices[x]] < _surfaceVertexSize);
  }

  computeSurfaceVertexOneRings();

  _surfaceTetIDs.clear();
  for(unsigned int x = 0; x < _tets.size(); x++)
    if(isSurfaceTet(x))
      _surfaceTetIDs.push_back(x);

  loadMaterials();

  _x.resize(_unconstrainedSize * 3);
  _x.setZero();

  _vertexNumberOfCopies.resize(_unconstrainedSize, 1);
    
  if(simulate){
    int rank = _unconstrainedSize * 3;
    _F.resize(_tets.size());
    _R.resize(rank);
    _R.setZero();

    if(verbose){
      cout << " Initialize each tet...";
      flush(cout);
    }
    #if USING_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(unsigned int x = 0; x < _tets.size(); x++){
      _tets[x].init();
    }
    if(verbose){
      cout << " done. " << endl;
    }
  }
  if(dynamic)
    computeMasses();
}

//////////////////////////////////////////////////////////////////////
// Construct a vector of the triangle faces on the mesh surface
//////////////////////////////////////////////////////////////////////
void TET_MESH::computeSurfaceFaces()
{
  Real verbose = SIMPLE_PARSER::getBool("verbose", false);
  // try to read in a precomputed cache
  if(readSurfaceFaceCache()) 
    return;
  
  if(verbose){
    cout << " Computing surface faces..."; flush(cout);
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
  for (unsigned int x = 0; x < _tets.size(); x++)
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

  _surfaceFaces.clear();
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
      _surfaceFaces.push_back(faces[0]);
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
        _surfaceFaces.push_back(faces[x]);
        faceBelongsToTet.push_back(tetHash[iter->first][x]);
        whichFaceInTet.push_back(triangleHash[iter->first][x]);
      }
    }
    // if (verbose && i % (int)(hashSize / 10) == 0)
    // {
    //     cout << 100 * ((Real)i / hashSize) << "% ";
    //     flush(cout);
    // }
  }
  cout << "done. " << endl;

  _restSurfaceFaceAreas.clear();
  _restSurfaceVertexAreas.clear();

  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE& face = _surfaceFaces[x];

    Real triArea = face.area();
    _restSurfaceFaceAreas[&face] = triArea;

    for(int y = 0; y < 3; y++){
      if(_restSurfaceVertexAreas.find(face.vertices[y]) == _restSurfaceVertexAreas.end())
        _restSurfaceVertexAreas[face.vertices[y]] = triArea / 3.0;
      else
        _restSurfaceVertexAreas[face.vertices[y]] += triArea / 3.0;
    }
  }
  _surfaceVertices.clear();
  for(map<VEC3F*, Real>::iterator iter = _restSurfaceVertexAreas.begin(); iter != _restSurfaceVertexAreas.end(); iter++){
    _surfaceVertices.push_back(iter->first);
  }

  // write the surface faces to a file
  string filename = string(_filename);
  filename += string(".surfacefaces");
  FILE* file = fopen(filename.c_str(), "wb");
  int totalFaces = faceBelongsToTet.size();
  fwrite((void*)&totalFaces, sizeof(int), 1, file);
  for (int x = 0; x < totalFaces; x++)
  {
    int indices[3];
    indices[0] = _vertexID[_surfaceFaces[x].vertices[0]];
    indices[1] = _vertexID[_surfaceFaces[x].vertices[1]];
    indices[2] = _vertexID[_surfaceFaces[x].vertices[2]];
    fwrite((void*)indices, sizeof(int), 3, file);
  }
  fclose(file);
  
  if(verbose){
    cout << "# of surface faces " << _surfaceFaces.size() << endl;
    // cout << " # of surface vertices " << _surfaceVertices.size() << endl;
    assert(_surfaceVertices.size() == _surfaceVertexSize);
  }
}

