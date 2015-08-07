/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef TET_BVH_H
#define TET_BVH_H

#include <SETTINGS.h>
#include <deformCD/aabb.h>
#include <geometry/TET_MESH.h>
#include <dtgrid/SPARSE_SDF.h>
#include <geometry/VERTEX_BVH.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace::std;

class TET_BVH_NODE {
public:
  TET_BVH_NODE();
  TET_BVH_NODE(unsigned int);  
  TET_BVH_NODE(unsigned int *lst, unsigned int lst_num, aabb* tetBoxes, VEC3F* tetCenters);

  ~TET_BVH_NODE();

  aabb& refit(TET_MESH* tetMesh, bool useLowresTets);

  bool collide(TET_MESH* tetMesh, const VEC3F& node, VEC3F& restPenetratingPosition, bool useLowresTets);

  void collide(TET_MESH* tetMesh, VERTEX_BVH_NODE* other, vector<pair<int, VEC3F> >& collidingNodes, bool useLowresTets);

  friend class TET_BVH;

private:
  aabb _box;
  unsigned int _id;
  TET_BVH_NODE* _left;
  TET_BVH_NODE* _right;
};
class TET_BVH {
public:
  TET_BVH(TET_MESH* tetMesh, SPARSE_SDF* sparseSDF, bool useLowresTets);
  TET_BVH(TET_MESH* tetMesh, SPARSE_SDF* sparseSDF, const vector<int>& tetIDs, bool useLowresTets);
  ~TET_BVH();

  // void visulization(int level);
  // void collide(bvh_tree *);
  // void self_collide();
  void collide(VERTEX_BVH* other, vector<pair<int, VEC3F> >& collidingNodes);
  void collide(vector<VEC3F>& nodes, vector<pair<int, VEC3F> >& collidingNodes);
  // void color();
  void refit();
  // bool get_contact(int i, unsigned int &id1, unsigned int &id2);
private:
  void init(vector<int>& tetIDs);

private:
  TET_BVH_NODE* _root;
  TET_MESH* _tetMesh;
  bool _useLowresTets;
};

#endif
