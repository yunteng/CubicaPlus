#ifndef PARTITIONED_TET_MESH_CUBATURE_TRAINER_H
#define PARTITIONED_TET_MESH_CUBATURE_TRAINER_H

#include <util/SIMPLE_PARSER.h>
#include <util/NNLS.h>
#include <geometry/SUBSPACE_TET_MESH.h>

template<class MATERIAL_CACHE, class BONE>
class PARTITIONED_TET_MESH_CUBATURE_TRAINER
{
public:
  PARTITIONED_TET_MESH_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~PARTITIONED_TET_MESH_CUBATURE_TRAINER();

  int gatherTrainingData();

  void clearSamples() { _sampleInternalForces.clear(); };

  inline int numberOfSamples() { return _numberOfSamples; };
  inline int sampleDimension() { return _sampleDimension; };
  inline int totalCandidates() { return _totalCandidates; };
  inline void setCurrentPartition(int p) { _currentPartition = p;}

  inline const vector<VECTOR>& samples() { return _sampleInternalForces; };
  VECTOR getCandidateQuantitiy(int sampleID);

  inline int getTrueID(int sampleID) { return _partitionedTets[sampleID]; }

  inline int maxKeyPoints() 
  { 
    int keyPoints = SIMPLE_PARSER::getInt("max internal key tets", 1000);
    keyPoints = keyPoints < _totalCandidates ? keyPoints : _totalCandidates;

    return keyPoints;
  }
  inline Real errorTolerance()
  {
    return SIMPLE_PARSER::getFloat("internal force error tolerance", 0.05);
  }
  inline int maxIteration()
  {
    return SIMPLE_PARSER::getInt("internal force cubature max iteration", 100);
  }
  inline vector<int>& initialKeyPoints()
  {
    return _tetMesh->partitionedKeyTets(_currentPartition);
  }
  inline vector<Real>& initialWeights()
  {
    return _tetMesh->partitionedKeyWeights(_currentPartition);
  }
  inline string cubatureFilename()
  {
    return _cubatureFilename;
  }
  inline void setCubatureFilename(const string& name)
  {
    _cubatureFilename = name;
  }


private:
  SUBSPACE_TET_MESH* _tetMesh;
  MATERIAL_CACHE* _materialCache;
  RIGGER<BONE>* _rigger;

  int _currentPartition;
  vector<int> _partitionedTets;

  vector<VECTOR> _sampleInternalForces;
  vector<VECTOR> _sampleIndividualTetInternalForces;
  vector<Real> _inverseMagnitudes;

  string _cubatureFilename;

  int _numberOfSamples;
  int _sampleDimension;
  int _totalCandidates;
};
#include "PARTITIONED_TET_MESH_CUBATURE_TRAINER.inl"
#endif
