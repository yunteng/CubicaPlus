/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef TRIANGLE_H
#define TRIANGLE_H

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

using namespace std;

class TRIANGLE
{
public:
  TRIANGLE(VEC3F& v0, VEC3F& v1, VEC3F& v2, VEC3F color = VEC3F(1, 1, 1));
  TRIANGLE(const TRIANGLE& triangle);

  bool operator==(const TRIANGLE &RHS) const;

  void draw();
  void draw(const vector<VEC3F>& colors);

  VEC3F center();
  Real area() { 
    _area = 0.5 * (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]).norm();
    return _area;   
  };
  const VEC3F& normal() { 
    _normal = (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]);
    _normal.normalize();
    return _normal; 
  };
  bool baryCenter(const VEC3F& point, VEC3F& lambda);
  VEC3F interpPosition(VEC3F& lambda);

  void boundingBox(VEC3F& mins, VEC3F& maxs);
  // intersect with line segment
  bool intersects(const VEC3F& start, const VEC3F& end);
  VEC3F projection(const VEC3F& point);

  Real maxEdgeLength(){
    Real maxl = (*vertices[0] - *vertices[1]).norm();
    maxl = max(maxl, (*vertices[0] - *vertices[2]).norm());
    maxl = max(maxl, (*vertices[2] - *vertices[1]).norm());
    return maxl;
  }
  VEC3F* vertex(int i) { return vertices[i]; };

public:
  VEC3F* vertices[3];
private:
  VEC3F _color;
  VEC3F _normal;
  Real _area;
};
#endif
