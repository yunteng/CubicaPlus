/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef PENALTY_SCF_CUBATURE_TRAINER_H
#define PENALTY_SCF_CUBATURE_TRAINER_H

#include <util/SIMPLE_PARSER.h>
#include <util/NNLS.h>
#include <geometry/SUBSPACE_TET_MESH.h>

template<class BONE>
class PENALTY_SCF_CUBATURE_TRAINER
{
public:
  PENALTY_SCF_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~PENALTY_SCF_CUBATURE_TRAINER();

  // load the training candidates, compute their 
  // subspace collision responses
  int gatherTrainingData(const string& poseID);

  void clearSamples() { _samplePenaltyForces.clear(); };

  inline int numberOfSamples() { return _numberOfSamples; };
  inline int sampleDimension() { return _sampleDimension; };
  inline int totalCandidates() { return _totalCandidates; };
  inline void setPartitions(int left, int right) { _leftPartition = left; _rightPartition = right;}

  inline const vector<VECTOR>& samples() { return _samplePenaltyForces; };

  // get the concatenated subspace force 
  // for a collision pair
  VECTOR getCandidateQuantitiy(int sampleID);

  inline int getTrueID(int sampleID) { return sampleID; }

  inline int maxKeyPoints() 
  { 
    return _totalCandidates;
  }
  inline Real errorTolerance()
  {
    return SIMPLE_PARSER::getFloat("scf cubature error tolerance", 0.05);
  }
  inline int maxIteration()
  {
    return SIMPLE_PARSER::getInt("scf cubature max iteration", 100);
  }
  inline vector<int>& initialKeyPoints()
  {
    return _keyPoints;
  }
  inline vector<Real>& initialWeights()
  {
    return _keyWeights;
  }
  inline string cubatureFilename()
  {
    return _cubatureFilename;
  }
  inline void setCubatureFilename(const string& name)
  {
    _cubatureFilename = name;
  }
  // some training poses may have the same 
  // relateive transformation between leftPartition
  // and rightPartition, group them together and 
  // only select the pose that has the maximum 
  // collision points for training 
  void groupPoses(const vector<string>& snapshotIdx, vector<string>& trainingPoses);
  
  // compute the relative transformation betwwen 
  // leftPartition and rightPartition, transform 
  // the surface position of rightPartition to the 
  // local coordinate system of leftPartition
  VECTOR getTransformedColumn(string poseID);
  // write out the cubatures, the writeCubatures 
  // function in NNHTP_CUBATURE_GENERATOR
  // wouldn't work
  void writePairwiseCubatures(const vector<int>& keyPoints, const vector<Real>& keyWeights, FILE* file);

private:
  // find which neighboring partition this vertex is closest to
  int closestNeighborPartition(int partition, int vertexID);
  // read the collision vertices and their force 
  // responces for a single pose
  int loadCollisionPoints(const string& posePrefix);


private:
  SUBSPACE_TET_MESH* _tetMesh;

  RIGGER<BONE>* _rigger;

  // the current partition pari we are
  // training on
  int _leftPartition;
  int _rightPartition;

  // helper structure to store all the
  // information of a vertex-traingle pair
  struct SELF_COLLISION_PAIR{
    int vertexID;
    int surfaceID;
    int vertexPartition;
    int trianglePartition;
    vector<VECTOR> collisionResponses;
    vector<VECTOR> reducedCollisionResponses;
  };
  vector<SELF_COLLISION_PAIR> _individualCollisionResponses;

  // reduced collision responses
  vector<VECTOR> _samplePenaltyForces;

  // the inverse magnitude of _samplePenaltyForces
  vector<Real> _inverseMagnitudes;

  // the ids of the cubature points
  vector<int> _keyPoints;
  vector<Real> _keyWeights;

  string _cubatureFilename;

  // number of force snapshots for a
  // single frame
  int _numberOfSamples;
  // left partition rank + right partition rank
  int _sampleDimension;
  int _totalCandidates;
};
#include "PENALTY_SCF_CUBATURE_TRAINER.inl"
#endif
