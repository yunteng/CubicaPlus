/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <geometry/TET_MESH.h>
#include <set>
#include <util/IO.h>
#if USING_OPENMP
#include <omp.h>
#endif

void TET_MESH::writeLowresEmbeddedMesh(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");
  if(file == NULL)
  {
    cout << " Cannot write low resolution embedded mesh at " << filename << endl;
    return;
  }
  int totalVertices = _lowresRestPose.size();
  int dummy = 0;
  // write out vertex array sizes
  fwrite((void*)&_lowSurfaceVertexSize, sizeof(int), 1, file);
  fwrite((void*)&totalVertices, sizeof(int), 1, file);
  fwrite((void*)&dummy, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < _lowresRestPose.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = _lowresRestPose[x][0];
    nodeDouble[1] = _lowresRestPose[x][1];
    nodeDouble[2] = _lowresRestPose[x][2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // write out tet vertex lists
  int totalTets = _lowresTets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = _lowresVertexID[_lowresTets[x].vertices[0]];
    indices[1] = _lowresVertexID[_lowresTets[x].vertices[1]];
    indices[2] = _lowresVertexID[_lowresTets[x].vertices[2]];
    indices[3] = _lowresVertexID[_lowresTets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }

  fclose(file);
}
bool TET_MESH::readLowresEmbeddedMesh(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  if(file == NULL){
    cout << " cannot read low resolution embedded mesh " << filename << endl;
    return false;
  }
  // assume the embedded mesh has the same scaling as this
  Real scale = SIMPLE_PARSER::getFloat("mesh scaling", 1.0);
  int surfaceVertices = 0;
  fread((void*)&surfaceVertices, sizeof(int), 1, file);
  _lowSurfaceVertexSize = surfaceVertices;

  int totalVertices = 0;
  fread((void*)&totalVertices, sizeof(int), 1, file);
  int constrained = 0;
  fread((void*)&constrained, sizeof(int), 1, file);
  totalVertices += constrained;

  // read in vertex positions
  _lowresVertices.resize(totalVertices);

  _lowresRestPose.resize(totalVertices);

  for(int x = 0; x < totalVertices; x++){
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    _lowresVertices[x] = node;
    _lowresRestPose[x] = node;
  }
  // read in tet vertex lists
  int totalTets;
  fread((void*)&totalTets, sizeof(int), 1, file);

  _lowresTets.clear();
  for(int x = 0; x < totalTets; x++)
  {
    int indices[4];
    fread((void*)&indices[0], sizeof(int), 4, file);

    _lowresTets.push_back(TET(_lowresVertices[indices[0]],
                        _lowresVertices[indices[1]],
                        _lowresVertices[indices[2]],
                        _lowresVertices[indices[3]]));
  }

  fclose(file);

  cout << " Lowres Embedding: " << totalVertices << " vertices and " << totalTets << " tets" << endl;

  // vertex addres to id
  _lowresVertexID.clear();
  for(int x = 0; x < totalVertices; x++){
    _lowresVertexID[&_lowresVertices[x]] = x;
  }
  // set corresponding rest versions
  for(unsigned x = 0; x < _lowresTets.size(); x++){
    for(int y = 0; y < 4; y++){
      int id = _lowresVertexID[_lowresTets[x].vertices[y]];
      _lowresTets[x].setRest(y, &_lowresRestPose[id]);
    }
  }

  string surfaceFilename = filename + ".surfacefaces";
  file = fopen(surfaceFilename.c_str(), "rb");
  assert(file != NULL);

  int totalSurfaceFaces;
  fread((void*)&totalSurfaceFaces, sizeof(int), 1, file);
  _lowresSurfaceFaces.clear();
  // _restSurfaceVertexAreas.clear();

  for (int x = 0; x < totalSurfaceFaces; x++)
  {
    int indices[3];
    fread((void*)indices, sizeof(int), 3, file);
    TRIANGLE face(_lowresVertices[indices[0]], _lowresVertices[indices[1]], _lowresVertices[indices[2]]);

    _lowresSurfaceFaces.push_back(face);
  }
  fclose(file);

  cout << "read in " << totalSurfaceFaces << " low res surface faces" << endl;

  return true;
}
bool TET_MESH::readLowresEmbeddingMap(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  if(file == NULL){
    cout << " cannot read low res embedding map " << filename << endl;
    return false;
  }
  int size = 0;
  fread((void*)&size, sizeof(int), 1, file);
  if(size != _lowresRestPose.size()){
    cout << " bad low res embedding map" << endl;
    fclose(file);
    return false;
  }
  _lowToHighID.resize(size);
  for(int x = 0; x < size; x++){
    int id = 0;
    fread((void*)&id, sizeof(int), 1, file);
    _lowToHighID[x] = id;
  }

  fclose(file);
  return true;
}
void TET_MESH::writeLowresEmbeddingMap(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");
  if(file == NULL){
    cout << " cannot write low res embedding map at " << filename << endl;
    return;
  }
  int size = _lowToHighID.size();
  fwrite((void*)&size, sizeof(int), 1, file);

  for(unsigned int x = 0; x < _lowToHighID.size(); x++){
    int id = _lowToHighID[x];
    fwrite((void*)&id, sizeof(int), 1, file);
  }

  fclose(file);
  return;
}
void TET_MESH::mapToLowresEmbedding()
{
  _lowToHighID.resize(_lowresRestPose.size());
  for(unsigned int x = 0; x < _lowresRestPose.size(); x++){
    Real minDist = 1e9;
    int matchIndex = -1;
    for(unsigned int y = 0; y < _restPose.size(); y++){
      Real dist = (_restPose[y] - _lowresRestPose[x]).norm();
      if(dist < minDist){
        minDist = dist;
        matchIndex = y;
      }
    }
    _lowToHighID[x] = matchIndex;
  }
}
void TET_MESH::updateLowresEmbedding()
{
  static bool firstTime = true;
  if(firstTime){
    for(int x = 0; x < _lowresRestPose.size(); x++){
      int id = _lowToHighID[x];
      _lowresRestPose[x] = _restPose[_lowToHighID[x]];
    }
    firstTime = false;
  }
  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for(int x = 0; x < _lowresVertices.size(); x++){
    int id = _lowToHighID[x];
    _lowresVertices[x] = _vertices[_lowToHighID[x]];
  }
}

bool TET_MESH::readMeshFile(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  
  if (file == NULL){
    cout << "tetmesh file " << filename << " could not be found!" << endl;
    return false;
  }

  readMeshFile(file);

  fclose(file);

  if(SIMPLE_PARSER::getBool("verbose", false)){
    cout << " Done reading mesh file. " << endl;

    Real minX = _vertices[0][0];
    Real maxX = _vertices[0][0];
    Real minY = _vertices[0][1];
    Real maxY = _vertices[0][1];
    Real minZ = _vertices[0][2];
    Real maxZ = _vertices[0][2];
    for (unsigned int x = 0; x < _vertices.size(); x++)
    {
      if (_vertices[x][0] < minX) minX = _vertices[x][0];
      if (_vertices[x][0] > maxX) maxX = _vertices[x][0];
      if (_vertices[x][1] < minY) minY = _vertices[x][1];
      if (_vertices[x][1] > maxY) maxY = _vertices[x][1];
      if (_vertices[x][2] < minZ) minZ = _vertices[x][2];
      if (_vertices[x][2] > maxZ) maxZ = _vertices[x][2];
    }

    cout << "Bounding box: " << endl
         << "(" << minX << ", " << minY << ", " << minZ << ")" << endl
         << "(" << maxX << ", " << maxY << ", " << maxZ << ")" << endl;
  }
  return true;

}
void TET_MESH::readMeshFile(FILE* file)
{
  Real scale = SIMPLE_PARSER::getFloat("mesh scaling", 1.0);
  Real verbose = SIMPLE_PARSER::getBool("verbose", false);
  // read in vertex array sizes
  fread((void*)&_surfaceVertexSize, sizeof(int), 1, file);
  fread((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fread((void*)&_constrainedSize, sizeof(int), 1, file);

  if (verbose){
    cout << " " << _surfaceVertexSize << " surface vertices" << endl 
         << " " << _unconstrainedSize  << " unconstrained nodes" << endl
         << " " << _constrainedSize << " constrained nodes" << endl;
  }

  _vertices.resize(_unconstrainedSize + _constrainedSize);

  if (verbose)
    cout << " Vertices size: " << (_unconstrainedSize + _constrainedSize) * sizeof(VEC3F) / pow(2.0, 20.0) << " MB " << endl;

  // read in vertex positions
  for (int x = 0; x < _unconstrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    _vertices[x] = node;
  }
  for (int x = 0; x < _constrainedSize; x++)
  {
    VEC3F node;
#ifdef SINGLE_PRECISION
    double nodeDouble[3];
    fread((void*)nodeDouble, sizeof(double), 3, file);
    node[0] = (float)nodeDouble[0];
    node[1] = (float)nodeDouble[1];
    node[2] = (float)nodeDouble[2];
#else
    fread((void*)&(node), sizeof(Real), 3, file);
#endif
    node *= scale;

    _vertices[x + _unconstrainedSize] = node;
  }

  // read in tet vertex lists
  int totalTets;
  fread((void*)&totalTets, sizeof(int), 1, file);

  if (verbose)
  {
    cout << " " << totalTets << " total tets" << endl
         << " Tets size: " << totalTets * sizeof(TET) / pow(2.0, 20.0) << " MB " << endl;
  }
  _tets.clear();
  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    fread((void*)&indices[0], sizeof(int), 4, file);

    _tets.push_back(TET(_vertices[indices[0]],
                        _vertices[indices[1]],
                        _vertices[indices[2]],
                        _vertices[indices[3]]));
  }
}

void TET_MESH::writeMeshFile(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");

  if (file == NULL){
    cout << "tetmesh file " << filename << " could not be opened!" << endl;
    return;
  }
  if(SIMPLE_PARSER::getBool("verbose", false))
    cout << " Writing file " << filename << " ... ";
  writeMeshFile(file);
  fclose(file);
  
  if(SIMPLE_PARSER::getBool("verbose", false))
    cout << " done. " << endl;
}
void TET_MESH::writeMeshFile(FILE* file)
{
  set<VEC3F*> surfaceVertexHash;
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++)
    surfaceVertexHash.insert(_surfaceVertices[x]);
  
  map<VEC3F*, int> newID;
  vector<VEC3F*> newVertices;
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++){
    newVertices.push_back(_surfaceVertices[x]);
    newID[_surfaceVertices[x]] = newVertices.size() - 1;
  }
  for(unsigned int x = 0; x < _vertices.size(); x++)
    if (surfaceVertexHash.find(&_vertices[x]) == surfaceVertexHash.end())
    {
      newVertices.push_back(&_vertices[x]);
      newID[&_vertices[x]] = newVertices.size() - 1;
    }

  // write out vertex array sizes
  fwrite((void*)&_surfaceVertexSize, sizeof(int), 1, file);
  fwrite((void*)&_unconstrainedSize, sizeof(int), 1, file);
  fwrite((void*)&_constrainedSize, sizeof(int), 1, file);

  // write out vertex positions
  for (unsigned int x = 0; x < newVertices.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = (*newVertices[x])[0];
    nodeDouble[1] = (*newVertices[x])[1];
    nodeDouble[2] = (*newVertices[x])[2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // write out tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = newID[_tets[x].vertices[0]];
    indices[1] = newID[_tets[x].vertices[1]];
    indices[2] = newID[_tets[x].vertices[2]];
    indices[3] = newID[_tets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }
}
//////////////////////////////////////////////////////////////////////
// write constrained mesh
//////////////////////////////////////////////////////////////////////
void TET_MESH::writeNewConstraints(vector<VEC3F*>& newConstraints)
{
  // make a new hash of the constrained nodes
  map<VEC3F*, bool> constrainedHash;
  for (unsigned int x = 0; x < newConstraints.size(); x++)
    constrainedHash[newConstraints[x]] = true;

  // record the translation from the old to new vertexIDs
  map<VEC3F*, int> newID;

  // total unconstrained and constrained nodes
  int unconstrainedSize = 0;
  int constrainedSize = 0;

  // go through all the vertices, record the unconstrained ones first
  vector<VEC3F*> newVertices;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (constrainedHash.find(&_vertices[x]) == constrainedHash.end())
    {
      newVertices.push_back(&_vertices[x]);
      newID[&_vertices[x]] = newVertices.size() - 1;
      unconstrainedSize++;
    }

  // now record the constrained ones
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (constrainedHash.find(&_vertices[x]) != constrainedHash.end())
    {
      newVertices.push_back(&_vertices[x]);
      newID[&_vertices[x]] = newVertices.size() - 1;
      constrainedSize++;
    }

  string filename = _filename + ".constrained";
  // start dumping to a new file
  FILE* file = fopen(filename.c_str(), "wb");
  
  if (file == NULL){
    cout << "Filename " << filename << " could not be opened!" << endl;
    return;
  }

  cout << " Writing file " << filename << " ... ";

  fwrite((void*)&_surfaceVertexSize, sizeof(int), 1, file);
  // write out vertex array sizes
  fwrite((void*)&unconstrainedSize, sizeof(int), 1, file);
  cout << "unconstrainedSize " << unconstrainedSize << endl;
  fwrite((void*)&constrainedSize, sizeof(int), 1, file);
  cout << "constrained " << constrainedSize << endl;

  // write out vertex positions
  for (unsigned int x = 0; x < newVertices.size(); x++)
  {
    double nodeDouble[3];
    nodeDouble[0] = (*newVertices[x])[0];
    nodeDouble[1] = (*newVertices[x])[1];
    nodeDouble[2] = (*newVertices[x])[2];
    fwrite((void*)nodeDouble, sizeof(double), 3, file);
  }

  // read in tet vertex lists
  int totalTets = _tets.size();
  fwrite((void*)&totalTets, sizeof(int), 1, file);

  for (int x = 0; x < totalTets; x++)
  {
    // for each node in the tet
    int indices[4];
    indices[0] = newID[_tets[x].vertices[0]];
    indices[1] = newID[_tets[x].vertices[1]];
    indices[2] = newID[_tets[x].vertices[2]];
    indices[3] = newID[_tets[x].vertices[3]];
    fwrite((void*)&indices[0], sizeof(int), 4, file);
  }

  fclose(file);
  cout << " done. " << endl;

  string faceFile = string(filename) + string(".surfacefaces");
  faceFile = string("rm ") + faceFile;
  system(faceFile.c_str());

}
//////////////////////////////////////////////
// read TetGen mesh and convert to our format
//////////////////////////////////////////////
bool TET_MESH::readTetGen(string prefix)
{
  Real scale = SIMPLE_PARSER::getFloat("mesh scaling", 1.0);
  bool verbose = SIMPLE_PARSER::getBool("verbose", false);
  string nodeFile = prefix + ".node";
  string tetFile = prefix + ".ele";
  string surfaceFile = prefix + ".face";

  fstream filestr;
  filestr.open(nodeFile.c_str(), fstream::in);
  if(!filestr.good()){
    cout << "failed to read " << nodeFile << endl;
    return false;
  }

  int totalNodes, dim;
  int attriCnt, boundaryMarker;

  filestr >> totalNodes >> dim >> attriCnt >> boundaryMarker;
  
  if(dim != 3){
    cout << "vertex dimension is not 3, format not supported!!!";
    return false;
  }

  if(verbose){
    cout << "reading in " << totalNodes << " vertices" << endl;

    cout << " Vertices size: " << totalNodes * sizeof(VEC3F) / pow(2.0, 20.0) << " MB " << endl;
  }

  _unconstrainedSize = totalNodes;
  _constrainedSize = 0;

  for(int x = 0; x < totalNodes; x++){
    int index;
    VEC3F node;
    filestr >> index >> node[0] >> node[1] >> node[2];
    node *= scale;
    _vertices.push_back(node);
    // do nothing about the attributes or the boundary marker for now
    int dummy;
    for(int y = 0; y < attriCnt; y++){
      filestr >> dummy;
    }
    if(boundaryMarker != 0){
      filestr >> dummy;
    }
  }
  filestr.close();
  
  
  // read in tet vertex lists
  filestr.open(tetFile.c_str(), fstream::in);
  if(!filestr.good()){
    cout << "failed to read " << tetFile << endl;
    return false;
  }
  int totalTets, nodesPerTet, attibute;
  filestr >> totalTets >> nodesPerTet >> attibute;
  if(nodesPerTet != 4){
    cout << "tet size is not 4, format not supported!!!" << endl;
    return false;
  }

  if(verbose){
    cout << totalTets << " total tets" << endl;
    cout << " Tets size: " << totalTets * sizeof(TET) / pow(2.0, 20.0) << " MB " << endl;
  }

  for (int x = 0; x < totalTets; x++)
  {
    int index, indices[4];
    int dummy;
    filestr >> index >> indices[0] >> indices[1] >> indices[2] >> indices[3];
    if(attibute == 1)
      filestr >> dummy;
    
    _tets.push_back(TET(_vertices[indices[0]],
                        _vertices[indices[1]],
                        _vertices[indices[2]],
                        _vertices[indices[3]]));
  }
  filestr.close();


  // read in surface triangle vertex lists
  filestr.open(surfaceFile.c_str(), fstream::in);
  if(filestr.good()){
    int totalSurfaceFaces;
    filestr >> totalSurfaceFaces >> attibute;

    if(verbose){
      cout << endl << totalSurfaceFaces << " total surface faces" << endl;
      cout << " Surface faces size: " << totalSurfaceFaces * sizeof(TRIANGLE) / pow(2.0, 20.0) << " MB " << endl;
    }

    for (int x = 0; x < totalSurfaceFaces; x++)
    {
      int index, indices[3];
      int dummy;
      filestr >> index >> indices[0] >> indices[1] >> indices[2];
      if(attibute == 1)
        filestr >> dummy;
      
      _surfaceFaces.push_back(TRIANGLE(_vertices[indices[0]],
                          _vertices[indices[1]],
                          _vertices[indices[2]]));
    }
    filestr.close();
  }
  // otherwise compute the surface faces from scratch
  else{
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
    }
  }
  
  set<VEC3F*> surfaceVertexSet;
  for (unsigned int x = 0; x < _surfaceFaces.size(); x++)
  {
    TRIANGLE& face = _surfaceFaces[x];
    for(int y = 0; y < 3; y++)
      surfaceVertexSet.insert(face.vertices[y]);
  }
  _surfaceVertices.clear();
  std::copy(surfaceVertexSet.begin(), surfaceVertexSet.end(), std::back_inserter(_surfaceVertices));
  _surfaceVertexSize = _surfaceVertices.size();

  return true;
}
//////////////////////////////////////////////////////////////////////
// read in the surface faces
// compute surface vertices
// compute surface triangle area 
// and surface vertex area at rest pose
//////////////////////////////////////////////////////////////////////
bool TET_MESH::readSurfaceFaceCache()
{
  if (_filename.size() == 0) return false;

  // write the surface faces to a file
  string filename = string(_filename);
  filename += string(".surfacefaces");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL) return false;
  
  int totalSurfaceFaces;
  fread((void*)&totalSurfaceFaces, sizeof(int), 1, file);
  _surfaceFaces.clear();
  _restSurfaceVertexAreas.clear();

  for (int x = 0; x < totalSurfaceFaces; x++)
  {
    int indices[3];
    fread((void*)indices, sizeof(int), 3, file);
    TRIANGLE face(_vertices[indices[0]], _vertices[indices[1]], _vertices[indices[2]]);

    _surfaceFaces.push_back(face);

    Real triArea = face.area();

    for(int y = 0; y < 3; y++){
      if(_restSurfaceVertexAreas.find(face.vertices[y]) == _restSurfaceVertexAreas.end())
        _restSurfaceVertexAreas[face.vertices[y]] = triArea / 3.0;
      else
        _restSurfaceVertexAreas[face.vertices[y]] += triArea / 3.0;
    }
  }
  _restSurfaceFaceAreas.clear();
  for(unsigned int x = 0; x < _surfaceFaces.size(); x++){
    _restSurfaceFaceAreas[&_surfaceFaces[x]] = _surfaceFaces[x].area();
  }

  fclose(file);
  cout << " surface faces cache found, # of surface faces " << _surfaceFaces.size() << endl;

  _surfaceVertices.clear();
  for(map<VEC3F*, Real>::iterator iter = _restSurfaceVertexAreas.begin(); iter != _restSurfaceVertexAreas.end(); iter++){
    _surfaceVertices.push_back(iter->first);
    int vid = _vertexID[iter->first];
  }

  cout << " # of surface vertices " << _surfaceVertices.size() << endl;

  return true;
}
void TET_MESH::loadMaterials()
{
  _totalMaterials = SIMPLE_PARSER::getInt("total materials", 0);
  if (_totalMaterials == 0)
  {
    cout << " NO MATERIALS SPECIFIED!!!!" << endl;
    _materialCopies = NULL;
    return;
  }

  // read in the actual materials
  // allocate a thread-safe copy of material for each thread
  _materialCopies = new MATERIAL**[_totalCores];
  _materialCopies[0] = new MATERIAL*[_totalMaterials];
  for (int x = 0; x < _totalMaterials; x++)
  {
    // read in the config file name for the material
    string materialString("material " + IO::intToString(x));

    string materialFilename = SIMPLE_PARSER::getString(materialString, "");

    // get the material
    MATERIAL* material = SIMPLE_PARSER::readMaterial(materialFilename);
    _materialCopies[0][x] = material;
  }

  // allocate a thread-safe copy of material for each thread
  // _materialCopies = new MATERIAL**[_totalCores];
  for (int x = 1; x < _totalCores; x++)
  {
    _materialCopies[x] = new MATERIAL*[_totalMaterials];
    for (int y = 0; y < _totalMaterials; y++)
      _materialCopies[x][y] = _materialCopies[0][y]->copy();
  }
  
  if(_totalMaterials > 1)
    readMaterialsIndices();
}
void TET_MESH::readMaterialsIndices()
{
  // construct the materials filename
  string filename = string(_filename);
  filename += string(".materials");
  FILE* file = fopen(filename.c_str(), "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " WARNING: No material distribution file found." << endl; 
    return;
  }
  int materialIndex;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    fread((void*)&materialIndex, sizeof(int), 1, file);
    // assert(materialIndex >= 0 && materialIndex < _totalMaterials);
    _tets[x].materialIndex() = materialIndex;
  }
  fclose(file);
}

void TET_MESH::writeMaterialIndices()
{
  // construct the materials filename
  string filename = string(_filename);
  filename += string(".materials");
  FILE* file = fopen(filename.c_str(), "wb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " WARNING: Could not write out materials file!" << endl; 
    return;
  }
  
  int materialIndex;
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    materialIndex = _tets[x].materialIndex();
    fwrite((void*)&materialIndex, sizeof(int), 1, file);
  }
  fclose(file);
}
void TET_MESH::readDisplacementFromRest(const string& filename)
{
  if(IO::read(_x, filename)){
    updateFullMesh();
  }
}
void TET_MESH::writeDisplacementFromRest(const string& filename)
{
  IO::write(_x, filename);
}

void TET_MESH::writeSurfaceVertexPositions(const string& filename)
{
  VECTOR sv(_surfaceVertices.size() * 3);
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++)
    sv.segment<3>(x * 3) = *_surfaceVertices[x];
  IO::write(sv, filename);
}
void TET_MESH::getSurfaceVertexDisp(VECTOR& sv)
{
  sv.resize(_surfaceVertices.size() * 3);
  for(unsigned int x = 0; x < _surfaceVertices.size(); x++){
    int id = _vertexID[_surfaceVertices[x]];
    sv.segment<3>(x * 3) = _vertices[id] - _restPose[id];
  }
}
void TET_MESH::writeObj(const string& filename)
{
  computeSurfaceVertexNormal();
  FILE* file = fopen(filename.c_str(), "w");
  if(file == NULL){
    cout << "Cannot open " << filename << " to write!!" << endl;
    return;
  }
  for(int x = 0; x < _surfaceVertexSize; x++){
    fprintf(file, "v %f %f %f\n", _vertices[x][0], _vertices[x][1], _vertices[x][2]);
  }

  for(int x = 0; x < _surfaceVertexSize; x++){
    fprintf(file, "vn %f %f %f\n", _surfaceNormals[x * 3], _surfaceNormals[x * 3 + 1], _surfaceNormals[x * 3 + 2]);
  }
  for(int x = 0; x < _surfaceFaces.size(); x++){
    TRIANGLE& face = _surfaceFaces[x];
    int indices[] = {_vertexID[face.vertices[0]],
                     _vertexID[face.vertices[1]],
                     _vertexID[face.vertices[2]]};
    fprintf(file, "f %d %d %d\n", indices[0] + 1, indices[1] + 1, indices[2] + 1);
  }
  fclose(file);
}
