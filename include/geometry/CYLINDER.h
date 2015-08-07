/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef CYLINDER_H
#define CYLINDER_H

#include <geometry/SURFACE.h>

class CYLINDER : public SURFACE  
{
public:
  CYLINDER();
  CYLINDER(Real* center, Real radius, Real* caps);
  virtual ~CYLINDER();

  virtual bool inside(const VEC3F& point);
  // virtual Real potential() { return 0.0f; };
  virtual Real distance(const VEC3F& point);
  virtual SURFACE* copy();
  virtual VEC3F contactPoint(const VEC3F& point);

  virtual void draw();

  // draw caps on the ends
  virtual void drawCapsule();

  QUATERNION& cylinderRotation() { return _cylinderRotation; };
  VEC3F& cylinderTranslation()   { return _cylinderTranslation; };
  const Real& radius()          { return _radius; };
  void setRadius(Real radius) { _radius = radius; _radiusSq = _radius * _radius; };
  Real getLength()             { return _caps[1] - _caps[0]; };
  void setLength(Real length);
  Real* caps()   { return _caps; };
  
  // return the bounding box dimensions
  virtual void boundingBox(VEC3F& mins, VEC3F& maxs);

  // get the (ODE VERSION) of the cylinder end caps, i.e.
  // aligned along the y axis
  void endPoints(VEC3F& leftCap, VEC3F& rightCap);

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);
  virtual MATRIX3 springJacobian(const VEC3F& collisionPoint);
  // virtual MATRIX3 dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);
  virtual void verifySpringJacobian(VEC3F& collisionPoint);

private:
  // cylinder rotation
  QUATERNION _cylinderRotation;

  // cylinder translation
  VEC3F _cylinderTranslation;

  // x,z coordinates of center axis
  Real _center[2];
  Real _radius;
  // bottom and top y coordinate of caps
  Real _caps[2];
  Real _radiusSq;
};

#endif
