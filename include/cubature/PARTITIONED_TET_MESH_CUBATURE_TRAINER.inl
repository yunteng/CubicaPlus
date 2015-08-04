/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <util/MATRIX_UTIL.h>

template<class MATERIAL_CACHE, class BONE>
PARTITIONED_TET_MESH_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::PARTITIONED_TET_MESH_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger),
  _currentPartition(-1),
  _numberOfSamples(0),
  _sampleDimension(0)
{
  if(rigger == NULL)
  {
    cout << " No skeleton embedding!!! Abort..." << endl;
    exit(0);
  }
  _materialCache = new MATERIAL_CACHE(tetMesh);
}
template<class MATERIAL_CACHE, class BONE>
PARTITIONED_TET_MESH_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::~PARTITIONED_TET_MESH_CUBATURE_TRAINER()
{

}

template<class MATERIAL_CACHE, class BONE>
int PARTITIONED_TET_MESH_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::gatherTrainingData()
{
  if(_currentPartition == -1){
    cout << " Haven't set the current training partition yet!!!" << endl;
    _numberOfSamples = 0;
    return _numberOfSamples;
  }

  _partitionedTets = _tetMesh->partitionedTets(_currentPartition);
  _totalCandidates = _partitionedTets.size();

  string dataPath = SIMPLE_PARSER::getString("data path", "");
  vector<string> snapshotIdx = IO::getAllSnapshots(dataPath);
  if(snapshotIdx.size() == 0){
    cout << " Could not find any training data !!!" << endl;
    return 0;
  }

  string posePath = SIMPLE_PARSER::getString("pose path", "");
  string skeletonPrefix = SIMPLE_PARSER::getString("skeleton prefix", "");

  _sampleInternalForces.resize(snapshotIdx.size());
  _sampleIndividualTetInternalForces.resize(snapshotIdx.size());

  int index = 0;

  for(unsigned int x = 0; x < snapshotIdx.size(); x++){
    // vector<MATRIX3> transformMatrix;

    string skeletonFilename = posePath + skeletonPrefix + snapshotIdx[x] + ".skeleton";

    _rigger->skeleton()->loadFrame(skeletonFilename);

    bool fromRest = true;

    TIMING_BREAKDOWN::tic();
    _rigger->updateSkinning(fromRest);

    vector<MATRIX3>& transformMatrix = _rigger->skinningRotation();

    _tetMesh->readDisplacementFromRest(dataPath + snapshotIdx[x] + ".state");

    _materialCache->cacheDecompositions();


    VECTOR individualForces;
    _materialCache->computeIndividualTetInternalForces(_partitionedTets, individualForces);

    _sampleInternalForces[index].resize(_tetMesh->partitionRank(_currentPartition));
    _sampleInternalForces[index].setZero();

    for(unsigned int y = 0; y < _partitionedTets.size(); y++)
    {
      TET& tet = _tetMesh->tets()[_partitionedTets[y]];

      for(int z = 0; z < 4; z++)
      {
        int vertexID = _tetMesh->vertexID(tet.vertices[z]);
        if(_tetMesh->isConstrained(vertexID))
          continue;

        int partitionedVertexID = _tetMesh->partitionedVertexID(_currentPartition, vertexID);

        individualForces.segment<3>(y * 12 + z * 3) = transformMatrix[vertexID].transpose() * individualForces.segment<3>(y * 12 + z * 3);

        _sampleInternalForces[index] += _tetMesh->partitionVertexBasis(_currentPartition, partitionedVertexID).transpose() * individualForces.segment<3>(y * 12 + z * 3);
      }
    }

    _sampleIndividualTetInternalForces[index] = individualForces;

    index++;
  
  }
  _sampleInternalForces.resize(index);
  _sampleIndividualTetInternalForces.resize(index);

  _inverseMagnitudes.resize(_sampleInternalForces.size());
  for(unsigned int x = 0; x < _sampleInternalForces.size(); x++){
    _inverseMagnitudes[x] = 1.0 / _sampleInternalForces[x].norm();
    _sampleInternalForces[x] *= _inverseMagnitudes[x];
  }

  _numberOfSamples = index;

  _sampleDimension = _sampleInternalForces[0].size();

  cout << "number of samples " << _numberOfSamples << endl;

  return _numberOfSamples;
}

template<class MATERIAL_CACHE, class BONE>
VECTOR PARTITIONED_TET_MESH_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::getCandidateQuantitiy(int sampleID)
{
  VECTOR trainingColumn(_numberOfSamples * _sampleDimension);
  trainingColumn.setZero();

  int tetID = _partitionedTets[sampleID];
  TET& tet = _tetMesh->tets()[tetID];

  for(unsigned int x = 0; x < _sampleIndividualTetInternalForces.size(); x++)
  {
    VECTOR fullForceSample = _sampleIndividualTetInternalForces[x].segment(sampleID * 12, 12);
    fullForceSample *= _inverseMagnitudes[x];
    for(int y = 0; y < 4; y++)
    {
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      if(_tetMesh->isConstrained(vertexID))
        continue;

      int partitionedVertexID = _tetMesh->partitionedVertexID(_currentPartition, vertexID);

      trainingColumn.segment(x * _sampleDimension, _tetMesh->partitionRank(_currentPartition)) += _tetMesh->partitionVertexBasis(_currentPartition, partitionedVertexID).transpose() * fullForceSample.segment<3>(y * 3);
    }
  }
  return trainingColumn;
}
