#include <util/MATRIX_UTIL.h>

template<class BONE>
PENALTY_SCF_CUBATURE_TRAINER<BONE>::PENALTY_SCF_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger),
  _leftPartition(-1),
  _rightPartition(-1),
  _numberOfSamples(0),
  _sampleDimension(0)
{
  if(rigger == NULL)
  {
    cout << " No skeleton embedding!!! Abort..." << endl;
    exit(0);
  }
}
template<class BONE>
PENALTY_SCF_CUBATURE_TRAINER<BONE>::~PENALTY_SCF_CUBATURE_TRAINER()
{

}

template<class BONE>
VECTOR PENALTY_SCF_CUBATURE_TRAINER<BONE>::getTransformedColumn(string poseID)
{
  string skeletonFilename = SIMPLE_PARSER::getString("pose path", "") + SIMPLE_PARSER::getString("skeleton prefix", "") + poseID + ".skeleton";

  _rigger->skeleton()->loadFrame(skeletonFilename);
  _rigger->skeleton()->fixSkeletonStructure();

  _tetMesh->readDisplacementFromRest(SIMPLE_PARSER::getString("data path", "") + poseID + ".state");

  VEC3F relativeTranslation;
  QUATERNION relativeRotation;
  _rigger->skeleton()->computeRelativeTransform(_leftPartition, _rightPartition, relativeTranslation, relativeRotation);

  VECTOR result(_tetMesh->partitionedSurfaceVertexSize(_rightPartition) * 3);
  result.setZero();
  for(int x = 0; x < _tetMesh->partitionedSurfaceVertexSize(_rightPartition); x++){
    int vertexID = _tetMesh->partitionedVertices(_rightPartition)[x];
    VEC3F pos = *(_tetMesh->vertex(vertexID));
    result.segment<3>(x * 3) = relativeRotation._transformVector(pos) + relativeTranslation;
  }

  return result;
}
template<class BONE>
void PENALTY_SCF_CUBATURE_TRAINER<BONE>::groupPoses(const vector<string>& snapshotIdx, vector<string>& trainingPoses)
{
  vector<vector<string> > poseIdGroups;

  vector<bool> grouped(snapshotIdx.size(), false);

  string posePath = SIMPLE_PARSER::getString("pose path", "");
  string skeletonPrefix = SIMPLE_PARSER::getString("skeleton prefix", "");
  
  for(unsigned int x = 0; x < snapshotIdx.size(); x++){
    if(grouped[x]) continue;
    
    vector<string> poseGroup;
    poseGroup.push_back(snapshotIdx[x]);

    string skeletonFilename = posePath + skeletonPrefix + snapshotIdx[x] + ".skeleton";

    _rigger->skeleton()->loadFrame(skeletonFilename);
    _rigger->skeleton()->fixSkeletonStructure();

    VEC3F relativeTranslation;
    QUATERNION relativeRotation;
    _rigger->skeleton()->computeRelativeTransform(_leftPartition, _rightPartition, relativeTranslation, relativeRotation);
    
    for(unsigned int y = x + 1; y < snapshotIdx.size(); y++){
      if(grouped[y]) continue;

      skeletonFilename = posePath + skeletonPrefix + snapshotIdx[y] + ".skeleton";
      _rigger->skeleton()->loadFrame(skeletonFilename);
      _rigger->skeleton()->fixSkeletonStructure();

      VEC3F tmpRelativeTranslation;
      QUATERNION tmpRelativeRotation;
      _rigger->skeleton()->computeRelativeTransform(_leftPartition, _rightPartition, tmpRelativeTranslation, tmpRelativeRotation);
      // if these pose y has the same configuration for partition left and right as pose x, group them together
      QUATERNION diffRot = tmpRelativeRotation - relativeRotation;
      if((tmpRelativeTranslation - relativeTranslation).squaredNorm() + diffRot.x() * diffRot.x() + diffRot.y() * diffRot.y() + diffRot.z() * diffRot.z() + diffRot.w() * diffRot.w() < 1e-4){
        poseGroup.push_back(snapshotIdx[y]);
        grouped[y] = true;
      }
      
    }
    grouped[x] = true;
    poseIdGroups.push_back(poseGroup);
  }

  for(unsigned int x = 0; x < poseIdGroups.size(); x++){
    int maxCollisionPoints = -1;
    int bestPoseId = -1;
    for(unsigned int y = 0; y < poseIdGroups[x].size(); y++){

      int totalCandidates = loadCollisionPoints(poseIdGroups[x][y]);
      if(totalCandidates > maxCollisionPoints){
        maxCollisionPoints = totalCandidates;
        bestPoseId = y;
      }
    }
    trainingPoses.push_back(poseIdGroups[x][bestPoseId]);
  }
  cout << "turned " << snapshotIdx.size() << " samples into " << trainingPoses.size() << " groups" << endl;
}
// find which neighboring partition this vertex is closest to
template<class BONE>
int PENALTY_SCF_CUBATURE_TRAINER<BONE>::closestNeighborPartition(int partition, int vertexID)
{
  VEC3F vertex = _tetMesh->restPose()[vertexID];
  Real minDist = 1e9;
  int closestNeighbor = -1;
  for(int x = 0; x < _rigger->skeleton()->totalBones(); x++){
    // only check neighbors
    if(x == partition || !_tetMesh->isNeighbors(x, partition))
      continue;
    VEC3F jointPosition = _rigger->skeleton()->computeRestJointPosition(x, partition);
    Real dist = (vertex - jointPosition).norm();
    if(dist < minDist){
      minDist = dist;
      closestNeighbor = x;
    }
  }
  // cout << "closest neighbor of vertex " << vertexID << " in parition " << partition << ": " << closestNeighbor << endl;
  return closestNeighbor;
}

template<class BONE>
int PENALTY_SCF_CUBATURE_TRAINER<BONE>::loadCollisionPoints(const string& posePrefix)
{
  _individualCollisionResponses.clear();
  string filename = SIMPLE_PARSER::getString("data path", "") + posePrefix + ".collisionresponse";
  FILE* file = fopen(filename.c_str(), "rb");
  if(file == NULL){
    return 0;
  }
  int size = 0;
  fread((void*)&size, sizeof(int), 1, file);

  for(int x = 0; x < size; x++){
    int vertexPartition = 0;
    int trianglePartition = 0;
    int vertexID = 0;
    int surfaceID = 0;
    int n = 0;
    fread((void*)&vertexPartition, sizeof(int), 1, file);
    fread((void*)&trianglePartition, sizeof(int), 1, file);
    fread((void*)&vertexID, sizeof(int), 1, file);
    fread((void*)&surfaceID, sizeof(int), 1, file);
    fread((void*)&n, sizeof(int), 1, file);

    // cout << x << " " << vertexPartition << " " << trianglePartition << " " << vertexID << " " << surfaceID << " " << n << endl;
    SELF_COLLISION_PAIR collisionPair;
    collisionPair.vertexPartition = vertexPartition;
    collisionPair.trianglePartition = trianglePartition;
    collisionPair.vertexID = vertexID;
    collisionPair.surfaceID = surfaceID;
    collisionPair.collisionResponses.resize(n);
    for(int y = 0; y < n; y++){
      collisionPair.collisionResponses[y].resize(12);
      for(int z = 0; z < 12; z++){
        Real val = 0;
        fread((void*)&val, sizeof(Real), 1, file);
        collisionPair.collisionResponses[y][z] = val;
      }
    }
    collisionPair.reducedCollisionResponses.clear();

    if(_leftPartition == vertexPartition && _rightPartition == trianglePartition)
      _individualCollisionResponses.push_back(collisionPair);
    else if(_rightPartition == vertexPartition && _leftPartition == trianglePartition)
      _individualCollisionResponses.push_back(collisionPair); 

    else if(vertexPartition == trianglePartition){
      if(vertexPartition == _leftPartition){
        if(closestNeighborPartition(vertexPartition, vertexID) == _rightPartition)
          _individualCollisionResponses.push_back(collisionPair);
      }else if(vertexPartition == _rightPartition){
        if(closestNeighborPartition(vertexPartition, vertexID) == _leftPartition)
          _individualCollisionResponses.push_back(collisionPair);
      }
    }

  }
  fclose(file);
  return _individualCollisionResponses.size();
} 

template<class BONE>
int PENALTY_SCF_CUBATURE_TRAINER<BONE>::gatherTrainingData(const string& poseID)
{
  if(_leftPartition == -1 || _rightPartition == -1){
    cout << " Haven't set the current training partitions yet!!!" << endl;
    _numberOfSamples = 0;
    return _numberOfSamples;
  }
  // cout << "left partition " << _leftPartition << " right partition " << _rightPartition << endl;

  loadCollisionPoints(poseID);
  _totalCandidates = _individualCollisionResponses.size();

  if(_totalCandidates == 0){
    _numberOfSamples = 0;
    return _numberOfSamples;
  }

  _numberOfSamples = _individualCollisionResponses[0].collisionResponses.size();

  string skeletonFilename = SIMPLE_PARSER::getString("pose path", "") + SIMPLE_PARSER::getString("skeleton prefix", "") + poseID + ".skeleton";

  _rigger->skeleton()->loadFrame(skeletonFilename);
  _rigger->skeleton()->fixSkeletonStructure();
  bool fromRest = true;

  TIMING_BREAKDOWN::tic();
  _rigger->updateSkinning(fromRest);

  vector<MATRIX3>& transformMatrix = _rigger->skinningRotation();

  _samplePenaltyForces.resize(_numberOfSamples);

  int leftRank = _tetMesh->partitionRank(_leftPartition);
  int rightRank = _tetMesh->partitionRank(_rightPartition);

  for(unsigned int x = 0; x < _numberOfSamples; x++){

    _samplePenaltyForces[x].resize(leftRank + rightRank);
    _samplePenaltyForces[x].setZero();

    for(unsigned int y = 0; y < _individualCollisionResponses.size(); y++)
    {
      VECTOR& collisionResponse = _individualCollisionResponses[y].collisionResponses[x];

      int vertexID = _individualCollisionResponses[y].vertexID;
      int vertexPartition = _individualCollisionResponses[y].vertexPartition;
      int trianglePartition = _individualCollisionResponses[y].trianglePartition;

      int surfaceID = _individualCollisionResponses[y].surfaceID;
      TRIANGLE& face = _tetMesh->surfaceFaces()[surfaceID];
      int triVID[3] = {_tetMesh->vertexID(face.vertices[0]),
                       _tetMesh->vertexID(face.vertices[1]),
                       _tetMesh->vertexID(face.vertices[2])
                      };
      
      bool valid = true;

      int partitionedVertexID = _tetMesh->partitionedVertexID(vertexPartition, vertexID);
      if(partitionedVertexID == -1)
        valid = false;
      int partitionedTriVID[3];

      for(int z = 0; z < 3 && valid; z++){
        partitionedTriVID[z] = _tetMesh->partitionedVertexID(trianglePartition, triVID[z]);
        if(partitionedTriVID[z] == -1){
          valid = false;
        }
      }

      VECTOR reducedCollisionResponse(leftRank + rightRank);
      reducedCollisionResponse.setZero();

      if(!valid){
        _individualCollisionResponses[y].reducedCollisionResponses.push_back(reducedCollisionResponse);
        continue;
      }

      collisionResponse.head<3>() = transformMatrix[vertexID].transpose() * collisionResponse.head<3>();
      for(int z = 0; z < 3; z++){
        collisionResponse.segment<3>((z + 1) * 3) = transformMatrix[triVID[z]].transpose() * collisionResponse.segment<3>((z + 1) * 3);
      }
      // cout << __LINE__ << endl;

      VECTOR reducedVertexForce = _tetMesh->partitionVertexBasis(vertexPartition, partitionedVertexID).transpose() * collisionResponse.head<3>();
      // cout << __LINE__ << endl;

      VECTOR reducedTriVertexForces[3];
      for(int z = 0; z < 3; z++){
        reducedTriVertexForces[z] = _tetMesh->partitionVertexBasis(trianglePartition, partitionedTriVID[z]).transpose() * collisionResponse.segment<3>((z + 1) * 3);
      }

      

      if(_leftPartition == vertexPartition){
        // cout << __LINE__ << endl;
        _samplePenaltyForces[x].head(leftRank) += reducedVertexForce;
        reducedCollisionResponse.head(leftRank) += reducedVertexForce;
      }else{

        _samplePenaltyForces[x].tail(rightRank) += reducedVertexForce;
        reducedCollisionResponse.tail(rightRank) += reducedVertexForce;
        // cout << __LINE__ << endl;
      }
      if(_leftPartition == trianglePartition){
        // cout << __LINE__ << endl;
        for(int z = 0; z < 3; z++){
          _samplePenaltyForces[x].head(leftRank) += reducedTriVertexForces[z];
          reducedCollisionResponse.head(leftRank) += reducedTriVertexForces[z];
        }
      }else{
        // cout << __LINE__ << endl;
        for(int z = 0; z < 3; z++){
          _samplePenaltyForces[x].tail(rightRank) += reducedTriVertexForces[z];
          reducedCollisionResponse.tail(rightRank) += reducedTriVertexForces[z];
        }
      }
      _individualCollisionResponses[y].reducedCollisionResponses.push_back(reducedCollisionResponse);
    }
  }

  _inverseMagnitudes.resize(_samplePenaltyForces.size());
  for(unsigned int x = 0; x < _samplePenaltyForces.size(); x++){
    _inverseMagnitudes[x] = 1.0 / _samplePenaltyForces[x].norm();
    _samplePenaltyForces[x] *= _inverseMagnitudes[x];
  }

  _sampleDimension = _samplePenaltyForces[0].size();

  return _numberOfSamples;
}

template<class BONE>
VECTOR PENALTY_SCF_CUBATURE_TRAINER<BONE>::getCandidateQuantitiy(int sampleID)
{
  VECTOR trainingColumn(_numberOfSamples * _sampleDimension);
  trainingColumn.setZero();

  vector<VECTOR>& reducedCollisionResponses = _individualCollisionResponses[sampleID].reducedCollisionResponses;
  for(unsigned int x = 0; x < reducedCollisionResponses.size(); x++)
    trainingColumn.segment(x * _sampleDimension, _sampleDimension) = reducedCollisionResponses[x] * _inverseMagnitudes[x];

  return trainingColumn;
}
template<class BONE>
void PENALTY_SCF_CUBATURE_TRAINER<BONE>::writePairwiseCubatures(const vector<int>& keyPoints, const vector<Real>& keyWeights, FILE* file)
{
  int cubatureSize = keyPoints.size();
  fwrite((void*)&cubatureSize, sizeof(int), 1, file);
  
  for(int x = 0; x < keyPoints.size(); x++){
    SELF_COLLISION_PAIR& scPair = _individualCollisionResponses[keyPoints[x]];
    fwrite((void*)&(scPair.vertexPartition), sizeof(int), 1, file);
    fwrite((void*)&(scPair.trianglePartition), sizeof(int), 1, file);
    fwrite((void*)&(scPair.vertexID), sizeof(int), 1, file);
    fwrite((void*)&(scPair.surfaceID), sizeof(int), 1, file);

    double dweight = keyWeights[x];
    fwrite((void*)&dweight, sizeof(double), 1, file);
    // cout << "vertex partition " << scPair.vertexPartition << " face partition " << scPair.trianglePartition << " vid " << scPair.vertexID << " fid " << scPair.surfaceID << " weight " << dweight << endl;
  }
}
