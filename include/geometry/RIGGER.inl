/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <set>
#include <algorithm>
#include <util/SIMPLE_PARSER.h>
#include <Eigen/SparseLU>
#include <util/IO.h>
#include <util/MATRIX_UTIL.h>
#include <geometry/DUAL_QUATERNION.h>
#if USING_OPENMP
#include <omp.h>
#endif

template<class BONE>
RIGGER<BONE>::RIGGER(SKELETON<BONE>* skeleton, TET_MESH* tetMesh):
  _skeleton(skeleton),
  _tetMesh(tetMesh)
{
  _skinningMethod = SIMPLE_PARSER::getString("skinning", "dual quaternion");
}
template<class BONE>
void RIGGER<BONE>::writeBoneWeights(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "wb");
  if(file == NULL){
    cout << "cannot open skinning file " << filename << " to write!!!" << endl;
    return;
  }
  int size = _skinning.size();
  fwrite((void*)&size, sizeof(int), 1, file);
  for(unsigned int x = 0; x < _skinning.size(); x++){
    int bn = _skinning[x].size();
    fwrite((void*)&bn, sizeof(int), 1, file);
    for(unsigned int y = 0; y < _skinning[x].size(); y++){
      int boneID = _skinning[x][y].first;
      double w = _skinning[x][y].second;
      fwrite((void*)&boneID, sizeof(int), 1, file);
      fwrite((void*)&w, sizeof(double), 1, file);
    }
  }
  fclose(file);
}

template<class BONE>
bool RIGGER<BONE>::readBoneWeights(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  if(file == NULL){
    cout << "cannot open skinning file " << filename << " to read!!!" << endl;
    return false;
  }
  int size;
  fread((void*)&size, sizeof(int), 1, file);
  _skinning.resize(size);
  _maxWeightIndex.resize(size);

  for(int x = 0; x < size; x++){
    int bn = 0;
    fread((void*)&bn, sizeof(int), 1, file);
    _skinning[x].resize(bn);

    int maxIndex = -1;
    Real maxW = -1;

    for(int y = 0; y < bn; y++){
      int boneID = 0;
      double w = 0;
      fread((void*)&boneID, sizeof(int), 1, file);
      fread((void*)&w, sizeof(double), 1, file);
      _skinning[x][y] = make_pair(boneID, (Real)w);
      if(w > maxW){
        maxW = w;
        maxIndex = boneID;
      }
    }
    _maxWeightIndex[x] = maxIndex;
  }
  fclose(file);
  return true;
}
template<class BONE>
void RIGGER<BONE>::drawBoneSkinning(int boneID)
{
  if(boneID < 0 || boneID > _skeleton->totalBones())
    return;

  const vector<VEC3F>& boneColors = _skeleton->colors();
  vector<TRIANGLE>& surfaceFaces = _tetMesh->surfaceFaces();
  for(unsigned int x = 0; x < surfaceFaces.size(); x++){
    int vertexIDs[3];
    VEC3F* vertices[3];
    for(int y = 0; y < 3; y++){
      vertices[y] = surfaceFaces[x].vertices[y];
      vertexIDs[y] = _tetMesh->vertexID(vertices[y]);
    }

    VEC3F vertexColors[3];
    bool draw = false;

    for (int y = 0; y < 3; y++)
    {
      vertexColors[y].setZero();

      vector<pair<int, Real> >& weightVector = _skinning[vertexIDs[y]];
      for (unsigned int z = 0; z < weightVector.size(); z++)
      {
        vertexColors[y] += boneColors[weightVector[z].first] * weightVector[z].second;

        if(boneID == weightVector[z].first)
          draw = true;
      }
    }
    if(!draw)
      continue;

    VEC3F normal = (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]);
    normal.normalize();

    VEC3F& v0 = *vertices[0];
    VEC3F& v1 = *vertices[1];
    VEC3F& v2 = *vertices[2];

    glBegin(GL_TRIANGLES);
      glNormal3d(normal[0], normal[1], normal[2]);
      glColor4f(vertexColors[0][0], vertexColors[0][1], vertexColors[0][2], 1);
      glVertex3f(v0[0], v0[1], v0[2]);
      glColor4f(vertexColors[1][0], vertexColors[1][1], vertexColors[1][2], 1);
      glVertex3f(v1[0], v1[1], v1[2]);
      glColor4f(vertexColors[2][0], vertexColors[2][1], vertexColors[2][2], 1);
      glVertex3f(v2[0], v2[1], v2[2]);
    glEnd();
  }
}
template<class BONE>
void RIGGER<BONE>::drawBoneSkinning()
{
  const vector<VEC3F>& boneColors = _skeleton->colors();
  vector<TRIANGLE>& surfaceFaces = _tetMesh->surfaceFaces();
  for(unsigned int x = 0; x < surfaceFaces.size(); x++){
    int vertexIDs[3];
    VEC3F* vertices[3];
    for(int y = 0; y < 3; y++){
      vertices[y] = surfaceFaces[x].vertices[y];
      vertexIDs[y] = _tetMesh->vertexID(vertices[y]);
    }

    // get the weights
    // float weights[3];
    VEC3F vertexColors[3];

    for (int y = 0; y < 3; y++)
    {
      vertexColors[y].setZero();

      vector<pair<int, Real> >& weightVector = _skinning[vertexIDs[y]];
      for (unsigned int z = 0; z < weightVector.size(); z++)
      {
        vertexColors[y] += boneColors[weightVector[z].first] * weightVector[z].second;
      }
    }

    VEC3F normal = (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]);
    normal.normalize();

    VEC3F& v0 = *vertices[0];
    VEC3F& v1 = *vertices[1];
    VEC3F& v2 = *vertices[2];

    glBegin(GL_TRIANGLES);
      glNormal3d(normal[0], normal[1], normal[2]);
      glColor4f(vertexColors[0][0], vertexColors[0][1], vertexColors[0][2], 1);
      glVertex3f(v0[0], v0[1], v0[2]);
      glColor4f(vertexColors[1][0], vertexColors[1][1], vertexColors[1][2], 1);
      glVertex3f(v1[0], v1[1], v1[2]);
      glColor4f(vertexColors[2][0], vertexColors[2][1], vertexColors[2][2], 1);
      glVertex3f(v2[0], v2[1], v2[2]);
    glEnd();
  }
}

/*
partition the tets based on the current skinning weights, each tet is associated with 
the most influencing bone
*/
template<class BONE>
void RIGGER<BONE>::buildSkinningPartition(vector<int>& tetPartitions)
{
  vector<TET>& tets = _tetMesh->tets();
  tetPartitions.resize(tets.size(), -1);
  for(unsigned int x = 0; x < tets.size(); x++){
    TET& tet = tets[x];
    map<int, int> partitionVotes;

    for(int y = 0; y < 4; y++){
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      partitionVotes[_maxWeightIndex[vertexID]]++;
    }

    int winningPartition = 0;
    int mostVotes = -1;
    map<int,int>::iterator iter;
    for (iter = partitionVotes.begin(); iter != partitionVotes.end(); iter++)
      if (iter->second > mostVotes)
      {
        mostVotes = iter->second;
        winningPartition = iter->first;
      } 

    tetPartitions[x] = winningPartition;
  }
}

/*
associate each vertex with its nearest bone
*/
template<class BONE>
void RIGGER<BONE>::buildRigidSkinning()
{
  cout << " Build rigid bone skinning ..."; flush(cout);
  // blow away any old skinning
  vector<VEC3F>& restPose = _tetMesh->restPose();
  _skinning.resize(restPose.size());

  _nearestBones.resize(restPose.size());
  _distsToNearest.resize(restPose.size());
  _maxWeightIndex.resize(restPose.size());

  vector<pair<VEC3F, VEC3F> > boneSegments;
  for(unsigned int x = 0; x < _skeleton->totalBones(); x++){
    boneSegments.push_back(_skeleton->bones()[x]->restBoneSegments());
  }

  // each vertex just picks up the skinning of the nearest surface vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    int closestBone = 0;
    Real closestDistance = 0;

    for (unsigned int y = 0; y < boneSegments.size(); y++)
    {
      VEC3F B = boneSegments[y].first;
      VEC3F P = restPose[x];
      VEC3F M = boneSegments[y].second - B;

      Real t0 = M.dot(P - B) / M.dot(M);

      Real distance = (P - (B + t0 * M)).norm();
      if (t0 <= 0)
        distance = (P - B).norm();

      if (t0 >= 1)
        distance = (P - (B + M)).norm();

      if (distance < closestDistance || y == 0)
      {
        closestBone = y;
        closestDistance = distance;
      }
    }
    _skinning[x].clear();
    _skinning[x].push_back(pair<int, Real>(closestBone, 1.0));

    _nearestBones[x] = closestBone;
    _maxWeightIndex[x] = closestBone;
    _distsToNearest[x] = closestDistance;
  }

  cout << "done." << endl;
}

/*
normalize the weight sum for each vertex to 1
*/
template<class BONE>
void RIGGER<BONE>::normalizeWeights()
{
  for(unsigned int x = 0; x < _skinning.size(); x++){
    Real sum = 0;
    if(_skinning[x].size() == 0)
      cout << "Vertex " << x << " has no weight set!!!" << endl;

    for(unsigned int y = 0; y < _skinning[x].size(); y++){
      sum += _skinning[x][y].second;
    }

    if(sum < 1e-7){
      cout << "Vertex " << x << " has weights sum to zero!!!" << endl;
    }else {
      Real scale = 1.0 / sum;
      for(unsigned int y = 0; y < _skinning[x].size(); y++){
        _skinning[x][y].second *= scale;
      }
    }
  }
}

/*
partition either the original mesh or its 
low-res embedding based on the current
skinning weights, used for collision detection
*/
template<class BONE>
void RIGGER<BONE>::buildSkinningPartition(vector<vector<int> >& partitionSurfaceVertices, vector<vector<int> >& partitionTets, bool useLowresTets)
{
  partitionSurfaceVertices.clear();
  partitionSurfaceVertices.resize(_skeleton->totalBones());

  if(_tetMesh->partitionedVertices().size() != _skeleton->totalBones()){
    vector<int> tetPartitions;
    buildSkinningPartition(tetPartitions);
    _tetMesh->buildPartitions(tetPartitions);
  }

  vector<vector<int> >& partitionedVertices = _tetMesh->partitionedVertices();

  for(unsigned int x = 0; x < partitionedVertices.size(); x++){
    for(int y = 0; y < partitionedVertices[x].size(); y++){
      if(_tetMesh->isSurfaceVertex(partitionedVertices[x][y]))
        partitionSurfaceVertices[x].push_back(partitionedVertices[x][y]);
    }
  }

  if(!useLowresTets){
    partitionTets = _tetMesh->partitionedTets();
    return;
  }

  vector<TET>& tets = _tetMesh->lowresTets();

  partitionTets.clear();
  partitionTets.resize(_skeleton->totalBones());
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    TET& tet = tets[x];
    map<int, int> partitionVotes;

    for (int y = 0; y < 4; y++)
    {
      // get the weight list of the current vertex
      int lowVertexID = _tetMesh->lowresVertexID(tet.vertices[y]);
      int vertexID = _tetMesh->lowToHighID()[lowVertexID];
      
      partitionVotes[_maxWeightIndex[vertexID]]++;
    }

    int winningPartition = 0;
    int mostVotes = -1;
    map<int,int>::iterator iter;
    for (iter = partitionVotes.begin(); iter != partitionVotes.end(); iter++)
      if (iter->second > mostVotes)
      {
        mostVotes = iter->second;
        winningPartition = iter->first;
      }
    partitionTets[winningPartition].push_back(x);
  }
}
/*
constrain tets that are penetraed by the skeleton
*/
template<class BONE>
void RIGGER<BONE>::constrainBoneTets()
{
  if(_tetMesh->constrained())
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Tet mesh is already constrained! " << endl;
    cout << " Don't try to constrain along bones twice, it will create redundant, confusing meshes! " << endl;
    return;
  }

  cout << " Constraining tets along bones ... "; flush(cout);

  vector<pair<VEC3F, VEC3F> > boneSegments;
  for(unsigned int x = 0; x < _skeleton->totalBones(); x++){
    boneSegments.push_back(_skeleton->bones()[x]->restBoneSegments());
  }

  // vertices to constrain
  map<VEC3F*, bool> constrainedVertices;

  // scan all tets, intersect with bones
  vector<TET>& tets = _tetMesh->tets();
  int constrainTets = 0;
  for (unsigned int x = 0; x < tets.size(); x++)
  {
    bool hitTet = false;

    for (int y = 0; y < _skeleton->totalBones() && !hitTet; y++)
    {
      // intersect line segment with each face
      for (int z = 0; z < 4; z++)
      {
        TRIANGLE face = tets[x].face(z);
        if(face.intersects(boneSegments[y].first, boneSegments[y].second))
        {
          hitTet = true;
          break;
        }
      }
      
    }

    // if the tet was hit, store its vertices
    if (hitTet){
      constrainTets++;
      for (int y = 0; y < 4; y++){
        if(!_tetMesh->isSurfaceVertex(tets[x].vertices[y]))
          constrainedVertices[tets[x].vertices[y]] = true;
      }
    }
  }
  cout << " done." << endl;
  cout << " Found " << constrainedVertices.size() << " vertices and " << constrainTets << " tets to constrain. " << endl;

  // flatten out to a vector
  vector<VEC3F*> newConstraints;
  for (map<VEC3F*, bool>::iterator iter = constrainedVertices.begin(); iter != constrainedVertices.end(); iter++)
    newConstraints.push_back(iter->first);

  cout << " Writing out new constrained mesh " << endl;
  _tetMesh->writeNewConstraints(newConstraints);
}

template<class BONE>
void RIGGER<BONE>::computeDQblend(int vertexID, DUAL_QUATERNION& dq_blend)
{
  vector<pair<int, Real> >& weights = _skinning[vertexID];

  QUATERNION maxq = _skeleton->bones()[_maxWeightIndex[vertexID]]->dualQuaternionFromRest().rotation();

  for (unsigned int y = 0; y < weights.size(); y++)
  {
    int b = weights[y].first;
    Real w = weights[y].second;

    DUAL_QUATERNION& dq = _skeleton->bones()[b]->dualQuaternionFromRest();

    if(dq.rotation().dot(maxq) < 0){
      w *= -1;
    }

    if(y == 0)
      dq_blend = dq * w;
    else
      dq_blend = dq_blend + (dq * w);
  }
  dq_blend.normalize();
}

template<class BONE>
void RIGGER<BONE>::updateDualQuaternionSkinning(bool fromRest)
{
  _skinningDisp.conservativeResize(_tetMesh->unconstrainedNodes() * 3);

  if(_skinningRotation.size() != _tetMesh->unconstrainedNodes())
    _skinningRotation.resize(_tetMesh->unconstrainedNodes());

  vector<VEC3F>& vertices = _tetMesh->vertices();
  vector<VEC3F>& restPose = _tetMesh->restPose();

  if(_previousDq_blend.size() != restPose.size())
    _previousDq_blend.resize(restPose.size());

  #if USING_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  // for each (constrained) vertex
  for (unsigned int x = 0; x < restPose.size(); x++)
  {
    DUAL_QUATERNION dq_blend;
    computeDQblend(x, dq_blend);

    if(fromRest)
      vertices[x] = dq_blend.transform(restPose[x]);
    else
      vertices[x] = dq_blend.transform(_previousDq_blend[x].inverseTransform(vertices[x]));

    _previousDq_blend[x] = dq_blend;

    if(x < _tetMesh->unconstrainedNodes()){
      _skinningDisp.segment<3>(x * 3) = dq_blend.transform(restPose[x]) - restPose[x];
      _skinningRotation[x] = dq_blend.rotation().toRotationMatrix();
    }
  }
}

/*
pull back the training samples to before-skinning space
*/
template<class BONE>
void RIGGER<BONE>::inverseTransformTrainingSamples()
{
  string dataPath = SIMPLE_PARSER::getString("data path", "");
  string posePath = SIMPLE_PARSER::getString("pose path", "");
  string skeletonPrefix = SIMPLE_PARSER::getString("skeleton prefix", "");

  vector<string> snapshotIdx = IO::getAllSnapshots(dataPath);
  
  vector<VECTOR> output;

  vector<VEC3F>& restPose = _tetMesh->restPose();

  for(unsigned int x = 0; x < snapshotIdx.size(); x++){
    VECTOR worldDisp;
    IO::read(worldDisp, dataPath + snapshotIdx[x] + ".state");

    string skeletonFile = posePath + skeletonPrefix + snapshotIdx[x] + ".skeleton";
    _skeleton->loadFrame(skeletonFile);

    VECTOR restDisp(worldDisp.size());
    int nVertices = restDisp.size() / 3;

    for(unsigned int y = 0; y < nVertices; y++){
      DUAL_QUATERNION dq_blend;
      computeDQblend(y, dq_blend);

      VEC3F translation = dq_blend.translation();
      restDisp.segment<3>(y * 3) = dq_blend.rotation().conjugate()._transformVector(worldDisp.segment<3>(y * 3) + restPose[y] - translation) - restPose[y];
    }
    IO::write(restDisp, dataPath + snapshotIdx[x] + ".restspace.state");
    output.push_back(restDisp);
  }
  MATRIX outputMat;
  MATRIX_UTIL::vectorToMatrix(output, outputMat, true);
  
  IO::write(outputMat, dataPath + "restspacedisp.matrix");
}
