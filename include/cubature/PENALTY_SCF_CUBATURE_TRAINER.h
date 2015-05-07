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

  int gatherTrainingData(const string& poseID);
  // int gatherTrainingData(vector<string>& snapshotIdx);
  void clearSamples() { _samplePenaltyForces.clear(); };

  inline int numberOfSamples() { return _numberOfSamples; };
  inline int sampleDimension() { return _sampleDimension; };
  inline int totalCandidates() { return _totalCandidates; };
  inline void setPartitions(int left, int right) { _leftPartition = left; _rightPartition = right;}

  inline const vector<VECTOR>& samples() { return _samplePenaltyForces; };
  VECTOR getCandidateQuantitiy(int sampleID);

  inline int getTrueID(int sampleID) { return sampleID; }

  inline int maxKeyPoints() 
  { 
    return _totalCandidates;
    // int keyPoints = SIMPLE_PARSER::getInt("max internal key tets", 1000);
    // keyPoints = keyPoints < _totalCandidates ? keyPoints : _totalCandidates;

    // return keyPoints;
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
    // string _cubatureFilename = _tetMesh->filename();
    // if(SIMPLE_PARSER::getBool("transformed basis", false))
    // {
    //   _cubatureFilename += ".transformed";
    // }
    // _cubatureFilename += ".internalForceCubature";
    return _cubatureFilename;
  }
  inline void setCubatureFilename(const string& name)
  {
    _cubatureFilename = name;
  }

  void groupPoses(const vector<string>& snapshotIdx, vector<string>& trainingPoses);
  
  VECTOR getTransformedColumn(string poseID);

  void writePairwiseCubatures(const vector<int>& keyPoints, const vector<Real>& keyWeights, FILE* file);

private:
  vector<string> getSnapshotsWithCollisions(const string& dataPath);
  int closestNeighborPartition(int partition, int vertexID);
  int loadCollisionPoints(const string& posePrefix);


private:
  SUBSPACE_TET_MESH* _tetMesh;

  RIGGER<BONE>* _rigger;

  int _leftPartition;
  int _rightPartition;

  struct SELF_COLLISION_PAIR{
    int vertexID;
    int surfaceID;
    int vertexPartition;
    int trianglePartition;
    vector<VECTOR> collisionResponses;
    vector<VECTOR> reducedCollisionResponses;
  };
  vector<SELF_COLLISION_PAIR> _individualCollisionResponses;

  vector<VECTOR> _samplePenaltyForces;
  // vector<VECTOR> _sampleIndividualTetInternalForces;
  vector<Real> _inverseMagnitudes;

  vector<int> _keyPoints;
  vector<Real> _keyWeights;

  string _cubatureFilename;

  int _numberOfSamples;
  int _sampleDimension;
  int _totalCandidates;
};
#include "PENALTY_SCF_CUBATURE_TRAINER.inl"
#endif
