#ifndef PLANE_H
#define PLANE_H

#include <geometry/SURFACE.h>

class PLANE : public SURFACE  
{
public:
  PLANE(const VEC3F& point, const VEC3F& normal);
  virtual ~PLANE();

  VEC3F& point() { return _point; };
  const VEC3F& normal() { return _normal; };
  MATRIX3& rotation() { return _rotation; };
  Real& adhesionStiffness() { return _adhesionStiffness; };
  Real& frictionCoeff() { return _frictionCoeff; };
  virtual void draw();
  virtual bool inside(const VEC3F& point);
  // signed distance < 0, point inside, > 0 outside
  virtual Real distance(const VEC3F& point);
  virtual VEC3F contactPoint(const VEC3F& point);

  virtual VEC3F force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity);

  virtual MATRIX3 springJacobian(const VEC3F& collisionPoint)
  {
    return _springJacobian;
  }

  // debug function -- verify that the spring Jacobian is correct
  virtual void verifySpringJacobian(VEC3F& collisionPoint);

  virtual void boundingBox(VEC3F& mins, VEC3F& maxs)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }

private:
  VEC3F _point;
  VEC3F _normal;

  Real _adhesionStiffness;
  Real _frictionCoeff;

  // when drawing the plane, apply this rotation about _point
  MATRIX3 _rotation;

  VEC3F _axis;
  Real _angle;

  MATRIX3 _springJacobian;
  MATRIX3 _dampingJacobian;
};

#endif
