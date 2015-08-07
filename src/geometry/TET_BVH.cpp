/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <geometry/TET_BVH.h>
#include <deformCD/aap.h>

extern float middle_xyz(char xyz, const VEC3F &p1, const VEC3F &p2, const VEC3F &p3);


TET_BVH::TET_BVH(TET_MESH* tetMesh, SPARSE_SDF* sparseSDF, bool useLowresTets):
 _tetMesh(tetMesh),
 _useLowresTets(useLowresTets)
{
  vector<TET>& tets = useLowresTets ? _tetMesh->lowresTets() : _tetMesh->tets();

  vector<int> narrowBandedTets;
  for(unsigned int x = 0; x < tets.size(); x++){
    // if any of the 4 vertices of the tet is in the narrow-banded SDF
    // add this tet to the BVH
    bool inBand = false;
    for(int y = 0; y < 4; y++){
      MyLocator dummyLoc;
      const VEC3F& vertex = *(tets[x].vertices[y]);
      Real pos[3];
      pos[0] = vertex[0];
      pos[1] = vertex[1];
      pos[2] = vertex[2];
      if(sparseSDF->getLocator(pos, &dummyLoc)){
        inBand = true;
        break;
      }
    }
    if(inBand)
      narrowBandedTets.push_back(x);
  }

  init(narrowBandedTets);
}

TET_BVH::TET_BVH(TET_MESH* tetMesh, SPARSE_SDF* sparseSDF, const vector<int>& tetIDs, bool useLowresTets):
 _tetMesh(tetMesh),
 _useLowresTets(useLowresTets)
{
  vector<TET>& tets = useLowresTets ? _tetMesh->lowresTets() : _tetMesh->tets();

  vector<int> narrowBandedTets;
  for(unsigned int x = 0; x < tetIDs.size(); x++){
    // if any of the 4 vertices of the tet is in the narrow-banded SDF
    // add this tet to the BVH
    bool inBand = false;
    for(int y = 0; y < 4; y++){
      MyLocator dummyLoc;
      const VEC3F& vertex = *(tets[tetIDs[x]].vertices[y]);
      Real pos[3];
      pos[0] = vertex[0];
      pos[1] = vertex[1];
      pos[2] = vertex[2];
      if(sparseSDF->getLocator(pos, &dummyLoc)){
        inBand = true;
        break;
      }
    }
    if(inBand)
      narrowBandedTets.push_back(tetIDs[x]);
  }

  init(narrowBandedTets);
}
void TET_BVH::init(vector<int>& tetIDs)
{
  aabb total;

  vector<TET>& tets = _useLowresTets ? _tetMesh->lowresTets() : _tetMesh->tets();
  for(unsigned int x = 0; x < tetIDs.size(); x++){
    for(int y = 0; y < 4; y++){
      total += *(tets[tetIDs[x]].vertices[y]);
    }
  }

  unsigned int totalTets = tetIDs.size();

  aabb* tetBoxes = new aabb[tets.size()];
  // aabb* tetBoxes = new aabb[totalTets];
  VEC3F* tetCenters = new VEC3F[tets.size()];

  aap  pln(total);
  unsigned int* idx_buffer = new unsigned int[totalTets];
  unsigned int left_idx = 0, right_idx = totalTets;

  for(unsigned int x = 0; x < totalTets; x++){
    unsigned int tetID = tetIDs[x];//narrowBandedTets[x];
    tetCenters[tetID] = tets[tetID].center();
    if(pln.inside(tetCenters[tetID]))
      idx_buffer[left_idx++] = tetID;
    else
      idx_buffer[--right_idx] = tetID;

    for(int y = 0; y < 4; y++){
      tetBoxes[tetID] += *(tets[tetID].vertices[y]);
    }
  }


  _root = new TET_BVH_NODE();
  _root->_box = total;
  if(left_idx == 0 || left_idx == totalTets){
    int hf = totalTets / 2;
    if(hf > 0){
      _root->_left = new TET_BVH_NODE(idx_buffer, hf, tetBoxes, tetCenters);
      _root->_right = new TET_BVH_NODE(idx_buffer + hf, totalTets - hf, tetBoxes, tetCenters);
    }else{
      _root->_left = new TET_BVH_NODE(idx_buffer, totalTets, tetBoxes, tetCenters);
      _root->_right = NULL;
    }
  }else{
    _root->_left = new TET_BVH_NODE(idx_buffer, left_idx, tetBoxes, tetCenters);
    _root->_right = new TET_BVH_NODE(idx_buffer + left_idx, totalTets - left_idx, tetBoxes, tetCenters);
  }

  delete[] tetBoxes;
  delete[] tetCenters;
  delete[] idx_buffer;
}
TET_BVH::~TET_BVH()
{
  delete _root;
}


void TET_BVH::refit()
{
  _root->refit(_tetMesh, _useLowresTets);
}

void TET_BVH::collide(vector<VEC3F>& nodes, vector<pair<int, VEC3F> >& collidingNodes)
{
  for(int x = 0; x < nodes.size(); x++){
    VEC3F restPenetratingPosition;
    if(_root->collide(_tetMesh, nodes[x], restPenetratingPosition, _useLowresTets))
      collidingNodes.push_back(make_pair(x, restPenetratingPosition));
  }
}
void TET_BVH::collide(VERTEX_BVH* other, vector<pair<int, VEC3F> >& collidingNodes)
{
  _root->collide(_tetMesh, other->_root, collidingNodes, _useLowresTets);
  for(unsigned x = 0; x < collidingNodes.size(); x++){
    int id = collidingNodes[x].first;
    collidingNodes[x].first = other->_vertexPIDs[id];
  }
}
//#################################################################

TET_BVH_NODE::TET_BVH_NODE()
{
  _id = UINT_MAX;
  _left = _right = NULL;
}

TET_BVH_NODE::TET_BVH_NODE(unsigned int id)
{
  _left = _right = NULL;
  _id = id;
}

TET_BVH_NODE::TET_BVH_NODE(unsigned int *lst, unsigned int lst_num, aabb* tetBoxes, VEC3F* tetCenters)
{
  assert(lst_num > 0);
  _left = _right = NULL;
  _id = UINT_MAX;
  
  if(lst_num == 1) {
    _id = lst[0];
    _box = tetBoxes[_id];
  }
  else { // try to split them
    for (unsigned int t = 0; t < lst_num; t++) {
      int i = lst[t];
      _box += tetBoxes[i];
    }

    if (lst_num == 2) { // must split it!
      _left = new TET_BVH_NODE(lst[0]);
      _right = new TET_BVH_NODE(lst[1]);
    } else {
      aap pln(_box);
      unsigned int left_idx = 0, right_idx = lst_num - 1;

      for (unsigned int t = 0; t < lst_num; t++) {
        int i = lst[left_idx];
        if (pln.inside(tetCenters[i]))
          left_idx++;
        else {// swap it
          unsigned int tmp = i;
          lst[left_idx] = lst[right_idx];
          lst[right_idx--] = tmp;
        }
      }

      int hal = lst_num / 2;
      if (left_idx == 0 || left_idx == lst_num)
      {
        _left = new TET_BVH_NODE(lst, hal, tetBoxes, tetCenters);
        _right = new TET_BVH_NODE(lst + hal, lst_num - hal, tetBoxes, tetCenters);
      }
      else {
        _left = new TET_BVH_NODE(lst, left_idx, tetBoxes, tetCenters);
        _right = new TET_BVH_NODE(lst + left_idx, lst_num - left_idx, tetBoxes, tetCenters);
      }

    }
  }
}

TET_BVH_NODE::~TET_BVH_NODE()
{
  if (_left) delete _left;
  if (_right) delete _right;
}

aabb& TET_BVH_NODE::refit(TET_MESH* tetMesh, bool useLowresTets)
{
  if(_left == NULL){
    _box.empty();
    for(int y = 0; y < 4; y++){
      TET& tet = useLowresTets ? tetMesh->lowresTets()[_id] : tetMesh->tets()[_id];
      _box += *(tet.vertices[y]);
    }
  }else {
    _box = _left->refit(tetMesh, useLowresTets);
    _box += _right->refit(tetMesh, useLowresTets);
  }

  return _box;
}
bool TET_BVH_NODE::collide(TET_MESH* tetMesh, const VEC3F& node, VEC3F& restPenetratingPosition, bool useLowresTets)
{
  if(!_box.inside(node))
    return false;
  if(_left == NULL){
    QUATERNION lambda;
    // tetMesh->tets()[_id].baryCenter(node, lambda, false);
    TET& tet = useLowresTets ? tetMesh->lowresTets()[_id] : tetMesh->tets()[_id];

    if(tet.baryCenter(node, lambda)){
      restPenetratingPosition = tet.interpRestPosition(lambda);
      return true;
    }else{
      // cout << "lambda " << lambda[0] << " " << lambda[1] << " " << lambda[2] << " " << lambda[3] << endl;
      return false;
    }
  }
  if(!_left->collide(tetMesh, node, restPenetratingPosition, useLowresTets))
    return _right->collide(tetMesh, node, restPenetratingPosition, useLowresTets);
  else
    return true;
}
void TET_BVH_NODE::collide(TET_MESH* tetMesh, VERTEX_BVH_NODE* other, vector<pair<int, VEC3F> >& collidingNodes, bool useLowresTets)
{
  if(!_box.overlaps(other->_box))
    return;

  if(_left == NULL && other->_left == NULL){
    VEC3F restPenetratingPosition;
    if(collide(tetMesh, other->_box._min, restPenetratingPosition, useLowresTets)){
      collidingNodes.push_back(make_pair(other->_id, restPenetratingPosition));
    }
    return;
  }
  if(_left == NULL){
    collide(tetMesh, other->_left, collidingNodes, useLowresTets);
    collide(tetMesh, other->_right, collidingNodes, useLowresTets);
  }else if(other->_left != NULL){
    _left->collide(tetMesh, other->_left, collidingNodes, useLowresTets);
    _left->collide(tetMesh, other->_right, collidingNodes, useLowresTets);
    _right->collide(tetMesh, other->_left, collidingNodes, useLowresTets);
    _right->collide(tetMesh, other->_right, collidingNodes, useLowresTets);
  }else{
    _left->collide(tetMesh, other, collidingNodes, useLowresTets);
    _right->collide(tetMesh, other, collidingNodes, useLowresTets);
  }
}
