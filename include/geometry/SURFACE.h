/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef SURFACE_H
#define SURFACE_H

#include <iostream>
#include <SETTINGS.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <util/SIMPLE_PARSER.h>

class SURFACE
{
public:
  SURFACE();
  
  virtual ~SURFACE();

  Real& collisionStiffness() { return _collisionStiffness; };
  Real& collisionDamping() { return _collisionDamping; };
  Real& thickness() { return _thickness; };

  virtual void draw() = 0;

  virtual bool inside(const VEC3F& point) = 0;
  virtual Real distance(const VEC3F& point) = 0;
  virtual VEC3F contactPoint(const VEC3F& point) = 0;

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity = VEC3F(0, 0, 0)) = 0;
  
  virtual MATRIX3 springJacobian(const VEC3F& collisionPoint) = 0;
  // debug function -- verify that the spring Jacobian is correct
  virtual void verifySpringJacobian(VEC3F& collisionPoint) = 0;

  // virtual MATRIX3 dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity = VEC3F(0, 0, 0));

  virtual void boundingBox(VEC3F& mins, VEC3F& maxs) = 0;

protected:
  string _type;
  Real _thickness;
  // collision response constants
  Real _collisionStiffness;
  Real _collisionDamping;
};

#endif
