#include <geometry/PLANE.h>
#include <util/MATRIX_UTIL.h>

PLANE::PLANE(const VEC3F& point, const VEC3F& normal) :
  _point(point), _normal(normal)
{
  _type.assign("PLANE");

  _normal.normalize();

  VEC3F xAxis(1,0,0);
  VEC3F zNew = xAxis.cross(_normal);
  zNew.normalize();

  // check for the case where the normal is the x axis
  if (fabs(fabs(xAxis.dot(_normal)) - 1.0) <=  1e-4)
    zNew = VEC3F(0,0,1);

  VEC3F xNew = _normal.cross(zNew);

  _rotation.col(0) = xNew;
  _rotation.col(1) = _normal;
  _rotation.col(2) = zNew;

  QUATERNION quaternion(_rotation);
  MATRIX_UTIL::quaternionToAxisAngle(quaternion, _axis, _angle);
  
  _adhesionStiffness = 0;
  _frictionCoeff = 0;

  _springJacobian = _normal * _normal.transpose();
  _springJacobian *= 0.5;
  _springJacobian += 0.5 * MATRIX3::Identity();
  _dampingJacobian = _springJacobian; 
  
  _springJacobian *= _collisionStiffness;
  _dampingJacobian *= _collisionDamping;
  
}

PLANE::~PLANE()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
bool PLANE::inside(const VEC3F& point) 
{
  Real dot = _normal.dot(point - _point);
  if (dot < 0.0) return true;

  // for resting contact, make points stick to the plane
  // if (dot < _stickiness) return true;

  return false;
}
VEC3F PLANE::contactPoint(const VEC3F& point)
{
  return point - (_normal.dot(point - _point)) * _normal;
}

Real PLANE::distance(const VEC3F& point) {
  VEC3F diff = point - _point;
  diff[0] = point[0] - _point[0];
  diff[1] = point[1] - _point[1];
  diff[2] = point[2] - _point[2];

  return diff.dot(_normal);
}

void PLANE::draw()
{
  VEC3F v0(-1, 0, -1);
  VEC3F v1(-1, 0, 1);
  VEC3F v2(1, 0, 1);
  VEC3F v3(1, 0, -1);

  glPushMatrix();
    glTranslatef(_point[0], _point[1], _point[2]);

    glRotatef(_angle, _axis[0], _axis[1], _axis[2]);

    glBegin(GL_TRIANGLES);
      glNormal3f(0,1,0);
      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
      glVertex3f(v2[0], v2[1], v2[2]);

      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v2[0], v2[1], v2[2]);
      glVertex3f(v3[0], v3[1], v3[2]);

    glEnd();
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glLineWidth(3.0f);
    glBegin(GL_LINES);
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.1f, 0.0f);
    glEnd();
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// compute the force on a point
//////////////////////////////////////////////////////////////////////
VEC3F PLANE::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
{
  // compute distance to plane
  /*VEC3F diff = collisionPoint - _point;

  Real distance = _normal * diff;
  //if (distance > 0) return VEC3F();

  Real velocityComponent = _normal * collisionVelocity;

  VEC3F springForce = _collisionStiffness * (distance * _normal);
  VEC3F dampingForce = _collisionDamping * (velocityComponent * _normal);

  return springForce + dampingForce;*/

  VEC3F diff = collisionPoint - _point;
  Real distance = _normal.dot(diff);
  VEC3F springForce = _springJacobian * (distance * _normal);
  return springForce;
  // VEC3F dampingForce = _dampingJacobian * collisionVelocity;
  // return springForce + dampingForce;
}

void PLANE::verifySpringJacobian(VEC3F& collisionPoint)
{
  VEC3F velocity;
  velocity.setZero();
  VEC3F originalForce = force(collisionPoint, velocity);
  VEC3F originalPoint = collisionPoint;

  cout << " original force: " << originalForce << endl;

  Real delta = 1e-6;
  MATRIX3 finiteDiff;
  for (int x = 0; x < 3; x++)
  {
    VEC3F perturb = originalPoint;
    perturb[x] += delta;

    VEC3F newForce = force(perturb, velocity);
    VEC3F diff = newForce - originalForce;
    diff *= 1.0 / delta;
    finiteDiff.col(x) = diff;
  }

  cout << " finite diff: " << finiteDiff << endl;
  cout << " computed: " << _springJacobian << endl;
}
