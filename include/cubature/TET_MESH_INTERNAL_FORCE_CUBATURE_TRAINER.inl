#include <util/MATRIX_UTIL.h>

template<class MATERIAL_CACHE, class BONE>
TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER(SUBSPACE_TET_MESH* tetMesh, RIGGER<BONE>* rigger):
  _tetMesh(tetMesh),
  _rigger(rigger),
  _numberOfSamples(0),
  _sampleDimension(0)
{
  _materialCache = new MATERIAL_CACHE(tetMesh);
}
template<class MATERIAL_CACHE, class BONE>
TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::~TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER()
{

}

template<class MATERIAL_CACHE, class BONE>
int TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::gatherTrainingData()
{
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

  for(unsigned int x = 0; x < snapshotIdx.size(); x++){
    vector<MATRIX3> transformMatrix;
    if(_rigger != NULL){
      string skeletonFilename = posePath + skeletonPrefix + snapshotIdx[x] + ".skeleton";

      _rigger->skeleton()->loadFrame(skeletonFilename);
      _rigger->skeleton()->fixSkeletonStructure();
      bool fromRest = true;

      TIMING_BREAKDOWN::tic();
      _rigger->updateSkinning(fromRest);
    }
    _tetMesh->readDisplacementFromRest(dataPath + snapshotIdx[x] + ".state");

    _materialCache->cacheDecompositions();

    if(_tetMesh->isSkinningBasis()){
      assert(_rigger != NULL);

      vector<MATRIX3>& transformMatrix = _rigger->skinningRotation();

      VECTOR& transformedForce = _materialCache->computeInternalForce();
      for(unsigned int y = 0; y < transformMatrix.size(); y++){
        transformedForce.segment<3>(y * 3) = transformMatrix[y].transpose() * transformedForce.segment<3>(y * 3);
      }
      _sampleInternalForces[x] = _tetMesh->U().transpose() * transformedForce;


      VECTOR individualForces;
      _materialCache->computeIndividualTetInternalForces(individualForces);

      vector<TET>& tets = _tetMesh->tets();
      for(unsigned int t = 0; t < tets.size(); t++){
        for(int s = 0; s < 4; s++){
          int vertexID = _tetMesh->vertexID(tets[t].vertices[s]);
          if(vertexID < transformMatrix.size()){
            individualForces.segment<3>(t * 12 + s * 3) = transformMatrix[vertexID].transpose() * individualForces.segment<3>(t * 12 + s * 3);
          }
          // ignore transform the forces for the constrained vertices, their basis will be 0 anyway
        }
      }
      _sampleIndividualTetInternalForces[x] = individualForces;

    }else{
      _sampleInternalForces[x] = _tetMesh->U().transpose() * _materialCache->computeInternalForce();
      _materialCache->computeIndividualTetInternalForces(_sampleIndividualTetInternalForces[x]);
    }    
  }

  _inverseMagnitudes.resize(_sampleInternalForces.size());
  for(unsigned int x = 0; x < _sampleInternalForces.size(); x++){
    _inverseMagnitudes[x] = 1.0 / _sampleInternalForces[x].norm();
    _sampleInternalForces[x] *= _inverseMagnitudes[x];
  }

  _numberOfSamples = snapshotIdx.size();
  _sampleDimension = _sampleInternalForces[0].size();

  return _numberOfSamples;
}

template<class MATERIAL_CACHE, class BONE>
VECTOR TET_MESH_INTERNAL_FORCE_CUBATURE_TRAINER<MATERIAL_CACHE, BONE>::getCandidateQuantitiy(int tetID)
{
  VECTOR trainingColumn(_numberOfSamples * _sampleDimension);

  MATRIX tetU = _tetMesh->tetBasis(tetID);
  for(unsigned int x = 0; x < _sampleIndividualTetInternalForces.size(); x++)
  {
    VECTOR fullForceSample = _sampleIndividualTetInternalForces[x].segment(tetID * 12, 12);
    fullForceSample *= _inverseMagnitudes[x];
    trainingColumn.segment(x * _sampleDimension, _sampleDimension) = tetU.transpose() * fullForceSample;
  }
  return trainingColumn;
}

