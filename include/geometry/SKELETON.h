/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef SKELETON_H
#define SKELETON_H

#include <SETTINGS.h>
#include <iostream>
#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <util/RGB_HSV.h>

using namespace::std;

template<class BONE>
class SKELETON
{
public:
  SKELETON();
  SKELETON(const string& filename);
  ~SKELETON();

  bool loadFrame(const string& filename);

  void writeFrame(const string& filename);

  void drawBones();

  /*
  load the hierarchy of skeleton
  */
  void loadSturcture(const string& filname);

  /*
  if the skeleton hierarchy is provided, adjust the current skeleton configuration so that the relative position of the two end points of each joint matches the rest configuration
  */
  void fixSkeletonStructure();
  
  int totalBones()              { return _bones.size();};
  vector<BONE*>& bones()        { return _bones; };
  const vector<VEC3F>& colors() { return _colors; }

  /*
  check if the current configuration is the same as rest pose
  */
  bool isRestPose()
  {
    for(unsigned int x = 0; x < _bones.size(); x++)
      if(!_bones[x]->isRestPose())
        return false;
    return true;
  }

  void computeRelativeTransform(int x, int y, VEC3F& relativeTranslation, QUATERNION& relativeRotation);
  VEC3F computeRestJointPosition(int x, int y);

private:
  /*
  Visualize the skinning weights
  */
  void computeColors();
  
private:
  vector<BONE*> _bones;
  vector<VEC3F> _colors;
  vector<vector<pair<int, Real> > > _skinning;

  /*
  top down ordering of the bones,
  _boneHierarchy[0] is the root
  */
  vector<int> _boneHierarchy;
};

#include "SKELETON.inl"

#endif
