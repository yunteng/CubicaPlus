#ifndef TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER_H
#define TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER_H

#include <util/SIMPLE_PARSER.h>
#include <util/NNLS.h>
#include <geometry/SUBSPACE_TET_MESH.h>

template<class MATERIAL_CACHE, class BONE>
class TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER
{
public:
  TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger = NULL);
  ~TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER();

  int gatherTrainingData();

  // relase some memory
  void clearSamples() { _sampleInternalForces.clear(); };

  int numberOfSamples() { return _numberOfSamples; };
  int sampleDimension() { return _sampleDimension; };
  int totalCandidates() { return _tetMesh->tets().size(); };
  const vector<VECTOR>& samples() { return _sampleInternalForces; };
  VECTOR getCandidateQuantitiy(int tetID);

  inline int getTrueID(int tetID) { return tetID; };

  inline int maxKeyPoints() 
  { 
    int keyPoints = SIMPLE_PARSER::getInt("max internal key tets", 1000);
    keyPoints = keyPoints < _tetMesh->tets().size() ? keyPoints : _tetMesh->tets().size();

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
    return _tetMesh->keyTets();
  }
  inline vector<Real>& initialWeights()
  {
    return _tetMesh->keyWeights();
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

  vector<VECTOR> _sampleInternalForces;
  vector<VECTOR> _sampleIndividualTetInternalForces;
  vector<Real> _inverseMagnitudes;

  string _cubatureFilename;

  int _numberOfSamples;
  int _sampleDimension;
};
#include "TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER.inl"
#endif
