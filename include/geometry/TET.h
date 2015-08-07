/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef TET_H
#define TET_H

#include <SETTINGS.h>
#include <iostream>
#include <geometry/TRIANGLE.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

class TET
{
public:
  TET(VEC3F& v0, VEC3F& v1, VEC3F& v2, VEC3F& v3);
  TET(const TET& tet);
  ~TET();

  void init();

  void setRest(int x, VEC3F* addr) { _restVertices[x] = addr; };
  Real restVolume() { return _restVolume; };
  MATRIX& PFPu() { return _PFPu; };
  // area vectors
  const VEC3F* b() { return _b; };
  inline const MATRIX3x4& bMat() { return _bMat; };
  // materal inverse
  const MATRIX3 DmInv() { return _DmInv; };

  inline int& materialIndex() { return _materialIndex; };
	inline const int materialIndex() const { return _materialIndex; };

  inline VEC3F interpPosition(const QUATERNION& lambda)
	{
	  return lambda.x() * (*vertices[0]) + lambda.y() * (*vertices[1]) + lambda.z() * (*vertices[2]) + lambda.w() * (*vertices[3]);
	}

	inline VEC3F interpRestPosition(const QUATERNION& lambda)
	{
	  return lambda.x() * (*_restVertices[0]) + lambda.y() * (*_restVertices[1]) + lambda.z() * (*_restVertices[2]) + lambda.w() * (*_restVertices[3]);
	}

  TRIANGLE face(int x);
  VEC3F center();
  VEC3F restCenter();
  Real volume();

  MATRIX3 F();
  MATRIX3 G();
  void computePFPu();

  bool baryCenter(const VEC3F& vert, QUATERNION& lambda);
  bool restBaryCenter(const VEC3F& vert, QUATERNION& lambda);

  void drawFaces();

public:
  VEC3F* vertices[4];

private:
  VEC3F* _restVertices[4];

  MATRIX _PFPu;
  MATRIX3 _DmInv;

  VEC3F _b[4];
  MATRIX3x4 _bMat;

  Real _restVolume;

  int _materialIndex;
};
#endif
