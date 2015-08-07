/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef VERTEX_BVH_H
#define VERTEX_BVH_H

#include <SETTINGS.h>
#include <deformCD/aabb.h>
#include <geometry/TET_MESH.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace::std;

class VERTEX_BVH_NODE {
public:
  VERTEX_BVH_NODE();
  VERTEX_BVH_NODE(unsigned int);  
  VERTEX_BVH_NODE(unsigned int *lst, unsigned int lst_num, VEC3F* centers);

  ~VERTEX_BVH_NODE();

  aabb& refit(const vector<VEC3F>& vertices);

  friend class VERTEX_BVH;
  friend class TET_BVH_NODE;

private:
  aabb _box;
  unsigned int _id;
  VERTEX_BVH_NODE* _left;
  VERTEX_BVH_NODE* _right;
};
class VERTEX_BVH {
public:
  VERTEX_BVH(TET_MESH* tetMesh);
  VERTEX_BVH(TET_MESH* tetMesh, const vector<int>& vertexIDs);
  ~VERTEX_BVH();

  // void visulization(int level);
  // void collide(bvh_tree *);
  // void self_collide();
  void collide(vector<VEC3F>& nodes, vector<pair<int, VEC3F> >& collidingNodes);
  // void color();
  void refit();
  // bool get_contact(int i, unsigned int &id1, unsigned int &id2);
  friend class TET_BVH;

private:
  void init();

private:
  TET_MESH* _tetMesh;
  vector<int> _vertexPIDs;
  vector<VEC3F> _vertices;
  VERTEX_BVH_NODE* _root;
};

#endif
