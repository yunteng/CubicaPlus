#include <geometry/SPHERE.h>
SPHERE::SPHERE(Real radius, const VEC3F& center):
  SURFACE(),
  _radius(radius),
  _center(center)
{
  _type.assign("SPHERE");
}
SPHERE::~SPHERE()
{
  
}

void SPHERE::draw()
{
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushMatrix();
    glTranslatef(_center[0], _center[1], _center[2]);
    glutSolidSphere(_radius, 100, 100);
  glPopMatrix();
}

bool SPHERE::inside(const VEC3F& point)
{
  return (point - _center).squaredNorm() <= (_radius + _thickness) * (_radius + _thickness);
}
Real SPHERE::distance(const VEC3F& point)
{
  return (point - _center).norm() - _radius;
}
VEC3F SPHERE::contactPoint(const VEC3F& point)
{
  VEC3F normal = point - _center;
  normal.normalize();
  return _center + (normal * _radius);
}

VEC3F SPHERE::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  VEC3F direction = collisionPoint - _center;
  VEC3F normal = direction;
  normal.normalize();

  Real velocityDot = normal.dot(collisionVelocity);

  VEC3F springForce = _collisionStiffness * (direction - (normal * _radius));
  VEC3F dampingForce = _collisionDamping * normal * velocityDot;

  return springForce + dampingForce;
}
MATRIX3 SPHERE::springJacobian(const VEC3F& collisionPoint)
{
  VEC3F direction = collisionPoint - _center;
  Real dot = direction.dot(direction);
  Real sqrtDot = sqrt(dot);
  //Real invSqrtDot = 1.0 / sqrtDot;
  MATRIX3 final;
  final.setZero();

  Real same = 1.0 - _radius / sqrtDot;
  final(0,0) = final(1,1) = final(2,2) = same;

  VEC3F scaledDirection = pow(dot, -1.5) * _radius * direction;

  VEC3F column0 = direction[0] * scaledDirection;
  final.col(0) += column0;

  VEC3F column1 = direction[1] * scaledDirection;
  final.col(1) += column1;

  VEC3F column2 = direction[2] * scaledDirection;
  final.col(2) += column2;

  final *= 0.5;
  final += 0.5 * MATRIX3::Identity();

  return _collisionStiffness * final;
}
void SPHERE::verifySpringJacobian(VEC3F& collisionPoint)
{
  VEC3F originalForce = force(collisionPoint);
  VEC3F originalPoint = collisionPoint;

  cout << " original force: " << originalForce << endl;

  Real delta = 1e-6;

  MATRIX3 finiteDiff;
  for (int x = 0; x < 3; x++)
  {
    VEC3F perturb = originalPoint;
    perturb[x] += delta;

    VEC3F newForce = force(perturb);
    VEC3F diff = newForce - originalForce;
    diff *= 1.0 / delta;
    finiteDiff.col(x) = diff;
  }
  
  cout << " finite diff: " << finiteDiff << endl;
  cout << " computed: " << springJacobian(collisionPoint) << endl;
}

// MATRIX3 SPHERE::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
// {

// }

bool SPHERE::intersect(const SPHERE& rightSphere)
{
  Real distance = (_center - rightSphere.center()).squaredNorm();
  Real combinedRadii = _radius + rightSphere.radius();
  return distance <= (combinedRadii * combinedRadii);
}

void SPHERE::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  mins = _center - VEC3F(_radius, _radius, _radius);
  maxs = _center + VEC3F(_radius, _radius, _radius);
}
