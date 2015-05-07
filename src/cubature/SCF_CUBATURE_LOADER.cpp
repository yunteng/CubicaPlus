#include <cubature/SCF_CUBATURE_LOADER.h>

#ifdef USING_OPENMP
#include <omp.h>
#endif
void SCF_CUBATURE::read(FILE* file)
{
  fread((void*)&vertexPartition, sizeof(int), 1, file);
  fread((void*)&trianglePartition, sizeof(int), 1, file);
  fread((void*)&vertexID, sizeof(int), 1, file);
  fread((void*)&surfaceID, sizeof(int), 1, file);
  double dWeight;
  fread((void*)&dWeight, sizeof(double), 1, file);
  weight = dWeight;
}
PAIRWISE_SCF_CUBATURES::PAIRWISE_SCF_CUBATURES():
  _isNeighbor(false),
  _cachedFoundMatch(false),
  _safeDistance(1.0)
{ 
  _cachedTranslation = VEC3F(1e3, 1e3, 1e3);
}
PAIRWISE_SCF_CUBATURES::~PAIRWISE_SCF_CUBATURES(){

}
VECTOR SCF_CUBATURE_LOADER::getTransformedColumn(int rightPartition, const VEC3F& relativeTranslation, const QUATERNION& relativeRotation)
{
  VECTOR result(_tetMesh->partitionedSurfaceVertexSize(rightPartition) * 3);

  for(unsigned int x = 0; x < _tetMesh->partitionedSurfaceVertexSize(rightPartition); x++){
    int vertexID = _tetMesh->partitionedVertices(rightPartition)[x];
    VEC3F pos = *(_tetMesh->vertex(vertexID));

    result.segment<3>(x * 3) = relativeRotation._transformVector(pos) + relativeTranslation;
  }
  return result;
}

bool PAIRWISE_SCF_CUBATURES::blendCubatures(VECTOR& transformedPosition, vector<SCF_CUBATURE>& output, bool isCheb)
{
  priority_queue<pair<int, Real>, vector<pair<int, Real> >, orderByValueGreaterThan> distanceHeap;
  
  VECTOR currentCoords = _pcaBasis.transpose() * transformedPosition;

  Real maxDist = 0;
  Real currentCoordsNorm = 1 / currentCoords.norm();

  Real minDist = 1e9;
  for(unsigned int x = 0; x < _sampleCoordinates.size(); x++){

    Real dist = (currentCoords - _sampleCoordinates[x]).norm() * currentCoordsNorm;

    maxDist = dist > maxDist ? dist : maxDist;
    minDist = dist < minDist ? dist : minDist;

    if(dist < 1e-3){
      output.insert(output.end(), _cubatures[x].begin(), _cubatures[x].end());
      return true;
    }else if(dist < 0.5){
      distanceHeap.push(make_pair(x, dist));
    }
  }

  if(distanceHeap.empty()) 
    return false;

  if(distanceHeap.size() == 1){
    pair<int, Real> neighborIdx = distanceHeap.top();
    vector<SCF_CUBATURE>& neighborCubatures = _cubatures[neighborIdx.first];
    output.insert(output.end(), neighborCubatures.begin(), neighborCubatures.end());
    return true;
  }else{
    pair<int, Real> neighbors[2];
    Real dist[2];
    Real weights[2];
    for(int x = 0; x < 2; x++){
      neighbors[x] = distanceHeap.top();
      dist[x] = neighbors[x].second;
      distanceHeap.pop();
    }
    if(dist[1] > 4 * dist[0]){
      vector<SCF_CUBATURE>& neighborCubatures = _cubatures[neighbors[0].first];
      output.insert(output.end(), neighborCubatures.begin(), neighborCubatures.end());

      return true;
    }

    weights[0] = dist[1] * dist[1] / (dist[0] * dist[0] + dist[1] * dist[1]);

    weights[1] = 1 - weights[0];
    
    map<pair<int, int>, SCF_CUBATURE> uniquePoints;
    
    for(int x = 0; x < 2; x++){
      vector<SCF_CUBATURE>& neighborCubatures = _cubatures[neighbors[x].first];

      for(unsigned int y = 0; y < neighborCubatures.size(); y++){
        SCF_CUBATURE cub = neighborCubatures[y];
        cub.weight *= weights[x];
        pair<int, int> idPair(cub.vertexPartition, cub.vertexID);
        if(uniquePoints.find(idPair) == uniquePoints.end())
          uniquePoints[idPair] = cub;
        else
          uniquePoints[idPair].weight += cub.weight;
      }
    }
    for(map<pair<int, int>, SCF_CUBATURE>::iterator iter = uniquePoints.begin(); iter != uniquePoints.end(); iter++){
      output.push_back(iter->second);
    }
    return true;
  }
  return false;
}


SCF_CUBATURE_LOADER::SCF_CUBATURE_LOADER(TET_MESH* tetMesh):
_tetMesh(tetMesh)
{
  if(_tetMesh->filename().find("cheb")==std::string::npos){
    _isCheb = false;
  }else{
    _isCheb = true;
  }
}
SCF_CUBATURE_LOADER::~SCF_CUBATURE_LOADER()
{

}

bool SCF_CUBATURE_LOADER::loadAllCubatures(const string& dirName)
{
  vector<string> files;
  int status = IO::getdir(dirName, files);
  if(status != 0){
    cout << __FILE__ << " " << __LINE__  << endl;
    cout << "failed to read scf cubatures from " << dirName << "!!!!!!!!" << endl;
    return false;
  }
  for(unsigned int x = 0; x < files.size(); x++){
    vector<string> elems;
    IO::split(files[x], '.', elems);
    vector<string>::iterator it = find(elems.begin(), elems.end(), "scfcubature");
    if(it == elems.end())
      continue;

    string pxStr = *(it + 1);
    string pyStr = *(it + 2);
    
    int px = atoi(pxStr.c_str());
    int py = atoi(pyStr.c_str());
    
    pair<int, int> partitionPair(px, py);

    _pairwiseCubatures[partitionPair].read(dirName + files[x]);

    _pairwiseCubatures[partitionPair].isNeighbor() = _tetMesh->isNeighbors(px, py);

    _pairwiseCubatures[partitionPair].leftPartition() = px;

    _pairwiseCubatures[partitionPair].rightPartition() = py;   
  }
  return true;
}
void PAIRWISE_SCF_CUBATURES::read(string filename)
{
  FILE* file = fopen(filename.c_str(), "rb");
  if(file == NULL){
    cout << "cannot open file " << filename << " to read" << endl;
    return;
  }
  IO::read(_pcaBasis, file);

  MATRIX coordsMat;
  IO::read(coordsMat, file);

  for(int x = 0; x < coordsMat.cols(); x++){
    _sampleCoordinates.push_back(coordsMat.col(x));
    vector<SCF_CUBATURE> cubaturePoints;
    int cnt;
    fread((void*)&cnt, sizeof(int), 1, file);

    for(int i = 0; i < cnt; i++){
      SCF_CUBATURE cubature;
      cubature.read(file);
      cubaturePoints.push_back(cubature);
    }
    _cubatures.push_back(cubaturePoints);
  }

  fclose(file);
}
bool SCF_CUBATURE_LOADER::getCubatureSurfaceVertices(pair<int, int> partitionPair, const VEC3F& relativeTranslation, const QUATERNION& relativeRotation, vector<pair<int, Real> >& leftVertices, vector<pair<int, Real> >& rightVertices)
{
  bool swap = false;
  if(partitionPair.first > partitionPair.second){
    swap = true;
    int tmp = partitionPair.first;
    partitionPair.first = partitionPair.second;
    partitionPair.second  = tmp;

  }
  map<pair<int, int>, PAIRWISE_SCF_CUBATURES>::iterator iter = _pairwiseCubatures.find(partitionPair);
  if(iter != _pairwiseCubatures.end()){

    if(iter->second.cacheValid(relativeTranslation, relativeRotation)){
      if(!(iter->second.cachedFoundMatch()))
        return false;
      if(swap){
        leftVertices = iter->second.cachedRightVertices();
        rightVertices = iter->second.cachedLeftVertices();
      }else{
        leftVertices = iter->second.cachedLeftVertices();
        rightVertices = iter->second.cachedRightVertices();
      }
      return true;
    }

    VECTOR transformedPosition = getTransformedColumn(partitionPair.second, relativeTranslation, relativeRotation);

    vector<SCF_CUBATURE> cubatures;

    bool foundMatch = iter->second.blendCubatures(transformedPosition, cubatures, _isCheb);
    iter->second.cachedFoundMatch() = foundMatch;
    if(!cubatures.empty()){
      if(!swap){
        for(unsigned int x = 0; x < cubatures.size(); x++){
          if(cubatures[x].vertexPartition == partitionPair.first)
            leftVertices.push_back(make_pair(cubatures[x].vertexID, cubatures[x].weight));
          else
            rightVertices.push_back(make_pair(cubatures[x].vertexID, cubatures[x].weight));
        }
        iter->second.cachedLeftVertices() = leftVertices;
        iter->second.cachedRightVertices() = rightVertices;
      }else{
        for(unsigned int x = 0; x < cubatures.size(); x++){
          if(cubatures[x].vertexPartition == partitionPair.second)
            leftVertices.push_back(make_pair(cubatures[x].vertexID, cubatures[x].weight));
          else
            rightVertices.push_back(make_pair(cubatures[x].vertexID, cubatures[x].weight));
        }
        iter->second.cachedLeftVertices() = rightVertices;
        iter->second.cachedRightVertices() = leftVertices;
      }
      
    }
    return foundMatch;
  }
  return false;
}
