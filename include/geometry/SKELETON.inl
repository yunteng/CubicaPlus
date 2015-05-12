#include <queue>

template<class BONE>
SKELETON<BONE>::SKELETON()
{
}

template<class BONE>
SKELETON<BONE>::SKELETON(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "r");
  if(file == NULL){
    cout << "No skeleton file " << filename << " found!" << endl;
    return;
  }

  int index = 0;
  while(!feof(file))
  {
    BONE* bone = new BONE(file);
    _bones.push_back(bone);
  }
  fclose(file);
  computeColors();
}

template<class BONE>
SKELETON<BONE>::~SKELETON()
{
  for(unsigned int x = 0; x < _bones.size(); x++)
    if(_bones[x]){
      delete _bones[x];
      _bones[x] = NULL;
    }
}
template<class BONE>
void SKELETON<BONE>::computeColors()
{
  _colors.resize(_bones.size());

  HSV begin(0, 1, 1);
  HSV end(320, 1, 1);
  if(_bones.size() == 1){
    _colors[0] = begin.toRGB().toVEC3F();
    return;
  }
  for(int x = 0; x < _bones.size(); x++){
    Real t = x * 1.0 / (_bones.size() - 1);
    _colors[x] = begin.lerp(end, t).toRGB().toVEC3F();
  }
}
template<class BONE>
void SKELETON<BONE>::loadSturcture(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "r");
  if(file == NULL){
    cout << " No skeleton structure file " << filename << " found!" << endl;
    return;
  }
  _boneHierarchy.clear();
  map<int, int> childToParent;
  while(!feof(file)){
    int child = 0;
    int parent = 0;
    fscanf(file, "%i %i\n", &child, &parent);
    childToParent[child] = parent;
  }
  if(childToParent.size() != _bones.size()){
    cout << " Invalid skeleton structure, abort!" << endl;
    return;
  }

  vector<bool> inHierarchy(_bones.size(), false);
  queue<int> Q;
  for(int x = 0; x < _bones.size(); x++)
    Q.push(x);

  while(!Q.empty()){
    int child = Q.front();
    Q.pop();
    if(inHierarchy[child])
      continue;
    if(childToParent[child] == -1 || inHierarchy[childToParent[child]]){
      inHierarchy[child] = true;
      _boneHierarchy.push_back(child);
    }else{
      Q.push(childToParent[child]);
      Q.push(child);
    }
  }
  assert(_boneHierarchy.size() == _bones.size());
  for(map<int, int>::iterator iter = childToParent.begin(); iter != childToParent.end(); iter++){
    _bones[iter->first]->setParentID(iter->second);
    if(iter->second != -1)
      _bones[iter->first]->setParent(_bones[iter->second]);
  }
  // cout << "bone hierarchy" << endl;
  // for(unsigned int x = 0; x < _boneHierarchy.size(); x++){
  //   cout << _boneHierarchy[x] << endl;
  // }
}

template<class BONE>
bool SKELETON<BONE>::loadFrame(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "r");
  if(file == NULL){
    cout << "No skeleton file " << filename << " found!" << endl;
    return false;
  }

  int index = 0;
  while(!feof(file))
  {
    if(_bones.size() < index + 1){
      BONE* bone = new BONE(file);
      _bones.push_back(bone);
    }else{
      _bones[index++]->updateTransform(file);
    }
  }
  fclose(file);

  if(SIMPLE_PARSER::getBool("verbose", false))
    cout << "loaded skeleton frame " << filename << endl;

  return true;
}

template<class BONE>
void SKELETON<BONE>::writeFrame(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "w");
  if(file == NULL){
    cout << " Cannot write skeleton to file " << filename <<"!" << endl;
    return;
  }

  int index = 0;
  for(unsigned int x = 0; x < _bones.size(); x++)
  {
    _bones[x]->write(file);
  }

  fclose(file);

  if(SIMPLE_PARSER::getBool("verbose", false))
    cout << "write skeleton frame " << filename << endl;
}

template<class BONE>
void SKELETON<BONE>::drawBones()
{
  glDisable(GL_DEPTH_TEST);
  for(unsigned int x = 0; x < _bones.size(); x++)
    _bones[x]->draw();
  glEnable(GL_DEPTH_TEST);
}

template<class BONE>
void SKELETON<BONE>::fixSkeletonStructure()
{
  for(unsigned int x = 0; x < _boneHierarchy.size(); x++){
    _bones[_boneHierarchy[x]]->updateTransformAccordingToParent();
  }
}

template<class BONE>
void SKELETON<BONE>::computeRelativeTransform(int x, int y, VEC3F& relativeTranslation, QUATERNION& relativeRotation)
{
  _bones[x]->computeRelativeRT(_bones[y], relativeTranslation, relativeRotation);
}
template<class BONE>
VEC3F SKELETON<BONE>::computeRestJointPosition(int x, int y)
{
  pair<VEC3F, VEC3F> leftEnds = _bones[x]->restBoneSegments();
  pair<VEC3F, VEC3F> rightEnds = _bones[y]->restBoneSegments();

  Real minDist = 1e9;
  VEC3F jointPosition;
  Real dist = (leftEnds.first - rightEnds.first).norm();
  if(dist < minDist){
    minDist = dist;
    jointPosition = 0.5 * (leftEnds.first + rightEnds.first);
  }
  dist = (leftEnds.first - rightEnds.second).norm();
  if(dist < minDist){
    minDist = dist;
    jointPosition = 0.5 * (leftEnds.first + rightEnds.second);
  }
  dist = (leftEnds.second - rightEnds.first).norm();
  if(dist < minDist){
    minDist = dist;
    jointPosition = 0.5 * (leftEnds.second + rightEnds.first);
  }
  dist = (leftEnds.second - rightEnds.second).norm();
  if(dist < minDist){
    minDist = dist;
    jointPosition = 0.5 * (leftEnds.second + rightEnds.second);
  }
  
  return jointPosition;
}
