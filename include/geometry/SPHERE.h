#ifndef SPHERE_H
#define SPHERE_H

#include <geometry/SURFACE.h>

class SPHERE: public SURFACE
{
public:
  SPHERE(Real radius = 0.5, const VEC3F& center = VEC3F(0, 0, 0));
  virtual ~SPHERE();

  Real& radius() { return _radius; };
  VEC3F& center() { return _center; };

  Real radius() const { return _radius; };
  const VEC3F& center() const { return _center; };

  virtual void draw();

  virtual bool inside(const VEC3F& point);
  // signed distance < 0, point inside, > 0 outside
  virtual Real distance(const VEC3F& point);
  virtual VEC3F contactPoint(const VEC3F& point);

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity = VEC3F(0, 0, 0));
  virtual MATRIX3 springJacobian(const VEC3F& collisionPoint);
  // debug function -- verify that the spring Jacobian is correct
  virtual void verifySpringJacobian(VEC3F& collisionPoint); 

  // virtual MATRIX3 dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity = VEC3F(0, 0, 0));

  virtual bool intersect(const SPHERE& rightSphere);

  virtual void boundingBox(VEC3F& mins, VEC3F& maxs);

private:
  Real _radius;
  VEC3F _center;
};

#endif
