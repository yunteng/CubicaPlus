#include <set>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <util/IO.h>
#include <util/MATRIX_UTIL.h>

SUBSPACE_TET_MESH::SUBSPACE_TET_MESH():
  TET_MESH(),
  _totalPartitionRank(0)
{
  _usePartitionedBasis = SIMPLE_PARSER::getBool("partitioned basis", false);

  if(_usePartitionedBasis)
    _isSkinningBasis = true;
  else
    _isSkinningBasis = SIMPLE_PARSER::getBool("transformed basis", false);

  /*
  load basis and cubatures
  */
  _basisFilename = _filename;
  if(_usePartitionedBasis){
    _basisFilename += ".transformed.partitionedbasis.matrix";
    return;
  }
  else if(_isSkinningBasis)
    _basisFilename += ".transformed.basis.matrix";
  else
    _basisFilename += ".basis.matrix";

  if(IO::read(_UBasis, _basisFilename)){
    
    cout << "Basis rank: " << _UBasis.cols() << endl;
    _q.resize(_UBasis.cols());
    _q.setZero();

    _keyTets.clear();
    _keyTetWeights.clear();
    // _keyTetBases.clear();

    if(!readCubatures())
      cout << " Did not find any internal force cubatures" << endl;
    else
      cout << " # of internal force cubatures " << _keyTets.size() << endl;

  }else{
    cout << " No subspace basis found!!" << endl;
  }
}

SUBSPACE_TET_MESH::~SUBSPACE_TET_MESH()
{

}

void SUBSPACE_TET_MESH::loadPartitionBases()
{
  cout << "Loading partitioned basis" << endl;
  _partitionBases.resize(_totalPartitions);
  _totalPartitionRank = 0;
  _partitionRankStartIdx.clear();
  
  for(int x = 0; x < _totalPartitions; x++)
  {
    _partitionRankStartIdx.push_back(_totalPartitionRank);
    IO::read(_partitionBases[x], _filename + ".transformed.partition." + IO::intToString(x) + ".matrix");
    _totalPartitionRank += _partitionBases[x].cols();

    cout << " partition " << x << " rank " << _partitionBases[x].cols() << endl;
  }
  cout << "Total rank " << _totalPartitionRank << endl;
  _q.resize(_totalPartitionRank);
  _q.setZero();

  _partitionedKeyTets.resize(_totalPartitions);
  _partitionedKeyWeights.resize(_totalPartitions);

  computeInterfaceUs();
  computeInternalInterfaceUs();
  cout << "Done." << endl;
}
void SUBSPACE_TET_MESH::computePartitionBases(Real varianceCutoff)
{
  _partitionBases.clear();
  _totalPartitionRank = 0;

  if(_partitionedVertices.size() == 0){
    cout << " No valid partitioning provided!!!" << endl;
    return;
  }

  string dataPath = SIMPLE_PARSER::getString("data path", "");
  MATRIX restdispMat;
  if(!IO::read(restdispMat, dataPath + "restspacedisp.matrix")){
    cout << " Cannot find any training data!!!" << endl;
    return; 
  }

  for(unsigned int x = 0; x < _partitionedVertices.size(); x++)
  {
    cout << " Computing PCA basis for partition " << x << endl;
    vector<int>& mask = _partitionedVertices[x];
    cout << "  dofs: " << mask.size() << endl;

    MATRIX localDispMat(mask.size() * 3, restdispMat.cols());
    localDispMat.setZero();
    for(int y = 0; y < restdispMat.cols(); y++){
      for(int z = 0; z < mask.size(); z++){
        localDispMat(z * 3, y) = restdispMat(mask[z
          ] * 3, y);
        localDispMat(z * 3 + 1, y) = restdispMat(mask[z] * 3 + 1, y);
        localDispMat(z * 3 + 2, y) = restdispMat(mask[z] * 3 + 2, y);
      }
    }

    MATRIX eigenvectors;
    VECTOR eigenvalues;
    MATRIX_UTIL::pcaEigen(localDispMat, false, eigenvectors, eigenvalues);

    Real sumAll = eigenvalues.sum();
    Real sumKeep = 0;

    int rank = 0;

    Real cutoffThreshold = varianceCutoff * sumAll;

    for(int i = 0; i < eigenvalues.size(); i++){
      sumKeep += eigenvalues[i];
      rank++;
      if(sumKeep >= cutoffThreshold)
        break;
    }

    if(rank < eigenvectors.cols())
      eigenvectors.conservativeResize(eigenvectors.rows(), rank);

    MATRIX_UTIL::orthogonalize(eigenvectors);

    cout << " Basis rank " << eigenvectors.cols() << endl;
    cout << " " << sumKeep / sumAll * 100 << "% of variance is retained" << endl;

    MATRIX_UTIL::verifyTruncatedPCA(localDispMat, eigenvectors);


    IO::write(eigenvectors, _filename + ".transformed.partition." + IO::intToString(x) + ".matrix");

    _partitionBases.push_back(eigenvectors);

    _totalPartitionRank += rank;
  }
  cout << "Total partition rank " << _totalPartitionRank << endl;
}

void SUBSPACE_TET_MESH::computePCABasis(Real varianceCutoff)
{
  string dataPath = SIMPLE_PARSER::getString("data path", "");
  vector<string> snapshotIdx = IO::getAllSnapshots(dataPath);
  if(snapshotIdx.size() == 0){
    cout << " Could not find any training data for computing PCA basis!!!" << endl;
    return;
  }
  vector<VECTOR> samples;
  string postfix(".state");
  if(_isSkinningBasis)
    postfix = ".restspace.state";

  for(unsigned int x = 0; x < snapshotIdx.size(); x++){
    VECTOR displacments;
    if(!IO::read(displacments, dataPath + snapshotIdx[x] + postfix)){
      cout << "IO error!!! abort... " << endl;
      exit(0);
    }
    samples.push_back(displacments);
  }
  cout << "performing PCA on " << samples.size() << " samples" << endl;

  MATRIX sampleMat;
  MATRIX_UTIL::vectorToMatrix(samples, sampleMat, true);

  VECTOR eigenvalues;
  MATRIX_UTIL::pcaEigen(sampleMat, false, _UBasis, eigenvalues);

  Real sumAll = eigenvalues.sum();
  Real sumKeep = 0;

  int rank = 0;

  Real cutoffThreshold = varianceCutoff * sumAll;

  for(int i = 0; i < eigenvalues.size(); i++){
    sumKeep += eigenvalues[i];
    rank++;
    if(sumKeep >= cutoffThreshold)
      break;
  }

  if(rank < _UBasis.cols())
  {
    _UBasis.conservativeResize(_UBasis.rows(), rank);
  }
  MATRIX_UTIL::orthogonalize(_UBasis);

  cout << " Basis rank " << _UBasis.cols() << endl;
  cout << " " << sumKeep / sumAll * 100 << "% of variance is retained" << endl;

  MATRIX_UTIL::verifyTruncatedPCA(sampleMat, _UBasis);

  _basisFilename = _filename;

  if(_isSkinningBasis)
    _basisFilename += ".transformed.basis.matrix";
  else
    _basisFilename += ".basis.matrix";

  IO::write(_UBasis, _basisFilename);
}

bool SUBSPACE_TET_MESH::readCubatures()
{
  string filename = _basisFilename + ".internalForceCubature";

  FILE* file = fopen(filename.c_str(), "rb");
   
  if(file == NULL){
    cout << " cannot open " << filename << " to read!!" << endl;
    return false;
  }

  int size = 0;
  fread((void*)&size, sizeof(int), 1, file);

  for(int x = 0; x < size; x++){
    int id = 0;
    double weight = 0;
    fread((void*)&id, sizeof(int), 1, file);
    fread((void*)&weight, sizeof(double), 1, file);
    _keyTets.push_back(id);
    _keyTetWeights.push_back(weight);
    // cout << "id " << id << " weight " << weight << endl;
    assert(id < _tets.size());
  }
  fclose(file);

  for(unsigned int x = 0; x < _keyTets.size(); x++)
      _tets[_keyTets[x]].init();

  return true;
}

bool SUBSPACE_TET_MESH::loadPartitionCubatures()
{
  _totalPartitionCubatures = 0;

  _partitionedKeyTets.resize(_totalPartitions);
  _partitionedKeyWeights.resize(_totalPartitionRank);

  for(unsigned int x = 0; x < _totalPartitions; x++){
    string filename = _filename + ".transformed.partitionedbasis.partition." + IO::intToString(x) + ".cubature";

    FILE* file = fopen(filename.c_str(), "rb");
   
    if(file == NULL){
      cout << " cannot open " << filename << " to read!!" << endl;
      return false;
    }
    _partitionedKeyTets[x].clear();
    _partitionedKeyWeights[x].clear();

    int size = 0;
    fread((void*)&size, sizeof(int), 1, file);

    for(int y = 0; y < size; y++){
      int id = 0;
      double weight = 0;
      fread((void*)&id, sizeof(int), 1, file);
      fread((void*)&weight, sizeof(double), 1, file);
      _partitionedKeyTets[x].push_back(id);
      _partitionedKeyWeights[x].push_back(weight);
      // assert(id < _tets.size());
      // cout << "id " << id << " weight " << weight << endl;
    }
    fclose(file);
    // _partitionedKeyTets[x] = _partitionedTets[x];
    // _partitionedKeyWeights[x].resize(_partitionedKeyTets[x].size(), 1.0);

    cout << " # of key tets for partition " << x << ": " << _partitionedKeyTets[x].size() << endl;

    _totalPartitionCubatures += _partitionedKeyTets[x].size();
  }
  cout << "total cubatures " << _totalPartitionCubatures << endl;

  _keyTets.clear();
  _keyTetWeights.clear();

  int totalKeyTets = 0;
  _partitionKeyTetStartIdx.clear();

  for(unsigned int x = 0; x < _totalPartitions; x++){
    _partitionKeyTetStartIdx.push_back(totalKeyTets);

    for(unsigned int y = 0; y < _partitionedKeyTets[x].size(); y++){
      _keyTets.push_back(_partitionedKeyTets[x][y]);
      _keyTetWeights.push_back(_partitionedKeyWeights[x][y]);
    }
    totalKeyTets += _partitionedKeyTets[x].size();
  }

  for(unsigned int x = 0; x < _keyTets.size(); x++)
    _tets[_keyTets[x]].init();

  return true;  
}

void SUBSPACE_TET_MESH::updateSubspaceMesh()
{
  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
  cout << " shouldn't be calling this for now!" << endl;
  // _x = _UBasis * _q;
  // updateFullMesh();
}

void SUBSPACE_TET_MESH::generateKeyTetsF()
{
  if(_F.size() != _keyTets.size())
    _F.resize(_keyTets.size());

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int x = 0; x < _keyTets.size(); x++){
    _F[x] = _tets[_keyTets[x]].F();
  }
}

void SUBSPACE_TET_MESH::drawKeyTets()
{
  glDisable(GL_DEPTH_TEST);

  glColor3f(1.0, 0.0, 0.0);
  for(unsigned int x = 0; x < _keyTets.size(); x++){
    _tets[_keyTets[x]].drawFaces();
  }

  glEnable(GL_DEPTH_TEST);
}


void SUBSPACE_TET_MESH::computeInterfaceUs()
{
  _interfaceUs.clear();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = _interfaceSprings.begin(); iter != _interfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;
    vector<pair<int, int> >& springs = iter->second;

    _interfaceUs[iter->first].resize(5);

    MATRIX& leftU = _interfaceUs[iter->first][0];
    MATRIX& rightU = _interfaceUs[iter->first][1];

    leftU.resize(springs.size() * 3, _partitionBases[leftPartition].cols());
    rightU.resize(springs.size() * 3, _partitionBases[rightPartition].cols());

    leftU.setZero();
    rightU.setZero();

    for(unsigned int x = 0; x < springs.size(); x++){
      int leftVID = springs[x].first;
      int rightVID = springs[x].second;

      leftU.block(x * 3, 0, 3, leftU.cols()) = partitionVertexBasis(leftPartition, leftVID);

      rightU.block(x * 3, 0, 3, rightU.cols()) = partitionVertexBasis(rightPartition, rightVID);
    }
    _interfaceUs[iter->first][2] = leftU.transpose() * leftU;
    _interfaceUs[iter->first][3] = leftU.transpose() * rightU;
    _interfaceUs[iter->first][4] = rightU.transpose() * rightU;
  }
}

void SUBSPACE_TET_MESH::computeInternalInterfaceUs()
{
  _internalInterfaceUs.clear();

  for(map<pair<int, int>, vector<pair<int, int> > >::iterator iter = _internalInterfaceSprings.begin(); iter != _internalInterfaceSprings.end(); iter++){
    int leftPartition = iter->first.first;
    int rightPartition = iter->first.second;
    vector<pair<int, int> >& springs = iter->second;

    _internalInterfaceUs[iter->first].resize(5);

    MATRIX& leftU = _internalInterfaceUs[iter->first][0];
    MATRIX& rightU = _internalInterfaceUs[iter->first][1];

    leftU.resize(springs.size() * 3, _partitionBases[leftPartition].cols());
    rightU.resize(springs.size() * 3, _partitionBases[rightPartition].cols());

    leftU.setZero();
    rightU.setZero();

    for(unsigned int x = 0; x < springs.size(); x++){
      int leftVID = springs[x].first;
      int rightVID = springs[x].second;

      leftU.block(x * 3, 0, 3, leftU.cols()) = partitionVertexBasis(leftPartition, leftVID);

      rightU.block(x * 3, 0, 3, rightU.cols()) = partitionVertexBasis(rightPartition, rightVID);
    }
    _internalInterfaceUs[iter->first][2] = leftU.transpose() * leftU;
    _internalInterfaceUs[iter->first][3] = leftU.transpose() * rightU;
    _internalInterfaceUs[iter->first][4] = rightU.transpose() * rightU;
  }
}
