#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <geometry/VERTEX_BVH.h>
#include <deformCD/aap.h>

extern float middle_xyz(char xyz, const VEC3F &p1, const VEC3F &p2, const VEC3F &p3);


VERTEX_BVH::VERTEX_BVH(TET_MESH* tetMesh):
 _tetMesh(tetMesh)
{
  _vertexPIDs.clear();
  _vertices.clear();

  vector<VEC3F*>& surfaceVertices = _tetMesh->surfaceVertices();
  for(unsigned int x = 0; x < surfaceVertices.size(); x++){
    int vertexPID = _tetMesh->vertexID(surfaceVertices[x]);

    _vertices.push_back(*surfaceVertices[x]);
    _vertexPIDs.push_back(vertexPID);
  }

  init();
}

VERTEX_BVH::VERTEX_BVH(TET_MESH* tetMesh, const vector<int>& vertexIDs):
 _tetMesh(tetMesh)
{
  _vertexPIDs.clear();
  copy(vertexIDs.begin(), vertexIDs.end(), back_inserter(_vertexPIDs));

  _vertices.clear();

  for(unsigned int x = 0; x < _vertexPIDs.size(); x++){
    _vertices.push_back(_tetMesh->vertices()[_vertexPIDs[x]]);
  }

  init();
}
void VERTEX_BVH::init()
{
  aabb total;

  unsigned int totalVertices = _vertexPIDs.size();

  for(unsigned int x = 0; x < totalVertices; x++)
    total += _vertices[x];

  aap  pln(total);
  unsigned int* idx_buffer = new unsigned int[totalVertices];
  unsigned int left_idx = 0, right_idx = totalVertices;

  VEC3F* centers = &_vertices[0];//new VEC3F[totalVertices];
  for(unsigned int x = 0; x < totalVertices; x++){
    // centers[x] = _vertices[x];
    if(pln.inside(_vertices[x]))
      idx_buffer[left_idx++] = x;
    else
      idx_buffer[--right_idx] = x;
  }

  _root = new VERTEX_BVH_NODE();
  _root->_box = total;
  if(left_idx == 0 || left_idx == totalVertices){
    int hf = totalVertices / 2;
    if(hf > 0){
      _root->_left = new VERTEX_BVH_NODE(idx_buffer, hf, centers);
      _root->_right = new VERTEX_BVH_NODE(idx_buffer + hf, totalVertices - hf, centers);
    }else{
      _root->_left = new VERTEX_BVH_NODE(idx_buffer, totalVertices, centers);
      _root->_right = NULL;
    }
  }else{
    _root->_left = new VERTEX_BVH_NODE(idx_buffer, left_idx, centers);
    _root->_right = new VERTEX_BVH_NODE(idx_buffer + left_idx, totalVertices - left_idx, centers);
  }

  delete[] idx_buffer;
}
VERTEX_BVH::~VERTEX_BVH()
{
  delete _root;
}


void VERTEX_BVH::refit()
{
  // vector<VEC3F>& vertices = _tetMesh->vertices();
  // UNCONSTRAINED_SUBSPACE_TET_MESH* subMesh = (UNCONSTRAINED_SUBSPACE_TET_MESH*)_tetMesh; 
  // MATRIX3 rotation = subMesh->rotationQuaternion().toExplicitMatrix3x3();
  // VEC3F translation = subMesh->rigidTranslation();
  vector<VEC3F>& vertices = _tetMesh->vertices();

  for(unsigned int x = 0; x < _vertexPIDs.size(); x++)
  {
    _vertices[x] = vertices[_vertexPIDs[x]];
  }

  _root->refit(_vertices);
}


//#################################################################

VERTEX_BVH_NODE::VERTEX_BVH_NODE()
{
  _id = UINT_MAX;
  _left = _right = NULL;
}

VERTEX_BVH_NODE::VERTEX_BVH_NODE(unsigned int id)
{
  _left = _right = NULL;
  _id = id;
}

VERTEX_BVH_NODE::VERTEX_BVH_NODE(unsigned int *lst, unsigned int lst_num, VEC3F* centers)
{
  assert(lst_num > 0);
  _left = _right = NULL;
  _id = UINT_MAX;
  
  if(lst_num == 1) {
    _id = lst[0];
    _box += centers[_id];
  }
  else { // try to split them
    for (unsigned int t = 0; t < lst_num; t++) {
      int i = lst[t];
      _box += centers[i];
    }

    if (lst_num == 2) { // must split it!
      _left = new VERTEX_BVH_NODE(lst[0]);
      _right = new VERTEX_BVH_NODE(lst[1]);
    } else {
      aap pln(_box);
      unsigned int left_idx = 0, right_idx = lst_num - 1;

      for (unsigned int t = 0; t < lst_num; t++) {
        int i = lst[left_idx];
        if (pln.inside(centers[i]))
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
        _left = new VERTEX_BVH_NODE(lst, hal, centers);
        _right = new VERTEX_BVH_NODE(lst + hal, lst_num - hal, centers);
      }
      else {
        _left = new VERTEX_BVH_NODE(lst, left_idx, centers);
        _right = new VERTEX_BVH_NODE(lst + left_idx, lst_num - left_idx, centers);
      }

    }
  }
}

VERTEX_BVH_NODE::~VERTEX_BVH_NODE()
{
  if (_left) delete _left;
  if (_right) delete _right;
}

aabb& VERTEX_BVH_NODE::refit(const vector<VEC3F>& vertices)
{
  if(_left == NULL){
    _box.empty();
    _box += vertices[_id];
  }else {
    _box = _left->refit(vertices);
    _box += _right->refit(vertices);
  }

  return _box;
}
