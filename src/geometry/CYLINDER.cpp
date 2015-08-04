/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <geometry/CYLINDER.h>
#include <util/MATRIX_UTIL.h>

CYLINDER::CYLINDER():
  SURFACE()
{
  _type.assign("CYLINDER");

  _center[0] = 0;
  _center[1] = 0;

  _caps[0] = -0.3;
  _caps[1] = 0.3;

  _radius = 0.05f;

  _radiusSq = _radius * _radius;

  _cylinderTranslation = VEC3F(0,0,0);
  _cylinderRotation = QUATERNION(1,0,0,0);
}

CYLINDER::CYLINDER(Real* center, Real radius, Real* caps):
  SURFACE()
{
  _type.assign("CYLINDER");

  _center[0] = center[0];
  _center[1] = center[1];

  _caps[0] = caps[0];
  _caps[1] = caps[1];

  _radius = radius;
  _radiusSq = _radius * _radius;
  
  _cylinderTranslation = VEC3F(0,0,0);
  _cylinderRotation = QUATERNION(1,0,0,0);
}

CYLINDER::~CYLINDER()
{

}

//////////////////////////////////////////////////////////////////////
// deep copy
//////////////////////////////////////////////////////////////////////
SURFACE* CYLINDER::copy() {
  CYLINDER* final = new CYLINDER(_center, _radius, _caps);
  final->cylinderTranslation() = _cylinderTranslation;
  final->cylinderRotation() = _cylinderRotation;

  return final;
}

//////////////////////////////////////////////////////////////////////
// inside outside function
//////////////////////////////////////////////////////////////////////
bool CYLINDER::inside(const VEC3F& pointConst) {
  // transform into local frame
  VEC3F point = pointConst;
  point -= _cylinderTranslation;
  point = _cylinderRotation.toRotationMatrix().transpose() * point;

  if (point[1] > _caps[1] || point[1] < _caps[0]) return false;

  Real diff[2];
  diff[0] = point[0] - _center[0];
  diff[1] = point[2] - _center[1];

  Real magnitude = diff[0] * diff[0] + diff[1] * diff[1];

  return magnitude <= _radiusSq;
}
VEC3F CYLINDER::force(const VEC3F& collisionPoint, const VEC3F& collisionVelocity) {
  /*VEC3F direction = collisionPoint - _center;
  VEC3F normal = direction;
  normal.normalize();

  Real velocityDot = normal * collisionVelocity;

  VEC3F springForce = _collisionStiffness * (direction - (normal * _radius));
  VEC3F dampingForce = _collisionDamping * normal * velocityDot;

  return springForce + dampingForce;*/
  VEC3F center = _cylinderTranslation;
  center[0] += _center[0];
  center[2] += _center[1];

  VEC3F normal = _cylinderRotation.toRotationMatrix().col(1);
  VEC3F direction = (collisionPoint - center) - normal * (normal.dot(collisionPoint - center));
  Real d = direction.norm();
  VEC3F springForce = _collisionStiffness * (1.0 - _radius / d) * direction;

  return springForce;
}
VEC3F CYLINDER::contactPoint(const VEC3F& point)
{
  VEC3F center = _cylinderTranslation;
  center[0] += _center[0];
  center[2] += _center[1];

  VEC3F normal = _cylinderRotation.toRotationMatrix().col(1);
  VEC3F direction = (point - center) - normal * (normal.dot(point - center));
  Real d = direction.norm();
  return point + ((_radius / d - 1.0) * direction);
}
MATRIX3 CYLINDER::springJacobian(const VEC3F& collisionPoint) {
  VEC3F normal = _cylinderRotation.toRotationMatrix().col(1);

  VEC3F center = _cylinderTranslation;
  center[0] += _center[0];
  center[2] += _center[1];
  
  VEC3F direction = (collisionPoint - center) - normal * (normal.dot(collisionPoint - center));
  
  Real d = direction.norm();

  MATRIX3 M = MATRIX3::Identity() - normal * normal.transpose();

  VEC3F pDpX = M * direction;

  Real same = 1.0 - _radius / d;

  VEC3F scaledDirection = (pow(d, -3) * _radius) * direction;

  MATRIX3 final = same * M + scaledDirection * pDpX.transpose();

  return _collisionStiffness * final;
}
void CYLINDER::verifySpringJacobian(VEC3F& collisionPoint)
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
  cout << " computed: " << springJacobian(collisionPoint) << endl;
}
// MATRIX3 CYLINDER::dampingJacobian(const VEC3F& collisionPoint, const VEC3F& collisionVelocity)
// {
//   return MATRIX(3, 3);
// }
//////////////////////////////////////////////////////////////////////
// distance field function
//////////////////////////////////////////////////////////////////////
Real CYLINDER::distance(const VEC3F& pointConst) {
  VEC3F point = pointConst;
  // transform into local frame
  point -= _cylinderTranslation;
  point = _cylinderRotation.toRotationMatrix().transpose() * point;

  // distance from the cylinder center
  Real diff[2];
  diff[0] = point[0] - _center[0];
  diff[1] = point[2] - _center[1];
  Real distance = sqrt(diff[0] * diff[0] + diff[1] * diff[1]);

  if(point[1] <= _caps[1] && point[1] >= _caps[0]){
    distance = _radius - distance;

    Real bottom = _caps[0] - point[1];
    Real top    = point[1] - _caps[1];
    return max(distance, max(bottom, top));
  }

  // if it is outside and projects in the y direction 
  // onto the top or bottom cap
  if (distance < _radius)
  {
    if (point[1] > _caps[1])
      return point[1] - _caps[1];

    return _caps[0] - point[1];
  }

  // else it must be closest to a cap edge
  if (point[1] > _caps[1])
    diff[1] = point[1] - _caps[1];
  else
    diff[1] = _caps[0] - point[1];

  diff[0] = distance - _radius;

  return sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
}

//////////////////////////////////////////////////////////////////////
// draw the surface to OpenGL
//////////////////////////////////////////////////////////////////////
void CYLINDER::drawCapsule()
{
  Real length = fabs(_caps[0] - _caps[1]);
  Real radius = _radius;

  glPushMatrix();
    VEC3F axis;
    Real angle;
    MATRIX_UTIL::quaternionToAxisAngle(_cylinderRotation, axis, angle);
    // _cylinderRotation.axisAngle(axis, angle);
    glTranslatef(_cylinderTranslation[0], _cylinderTranslation[1], _cylinderTranslation[2]);
    glRotatef(angle, axis[0], axis[1], axis[2]);

    Real cylHalfHeight = length / 2.0;

    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i < 17; i++)
    {
      Real angle = i / 16.0 * 2.0 * M_PI;
      Real ca = cos(angle);
      Real sa = sin(angle);

      VEC3F normal(ca, sa, 0);
      normal = _cylinderRotation.toRotationMatrix() * normal;
      glNormal3f(normal[0], normal[1], normal[2]);
 
      glVertex3f(radius * ca, cylHalfHeight, radius * sa);
      glVertex3f(radius * ca, -cylHalfHeight, radius * sa);
    }

    glEnd();

    // draw the caps
    glTranslated(0, cylHalfHeight,0);
    glutSolidSphere(radius, 16, 12);
    glTranslated(0, -2.0 * cylHalfHeight,0);
    glutSolidSphere(radius, 16, 12);

    // draw a line towards one of the end caps
    glDisable(GL_DEPTH_TEST);
    glColor4f(0,10,0,1);
    glLineWidth(10);
    glBegin(GL_LINES);
      glVertex3f(0, 0, 0);
      glVertex3f(0, cylHalfHeight, 0);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// draw the surface to OpenGL
//////////////////////////////////////////////////////////////////////
void CYLINDER::draw()
{
  Real length = fabs(_caps[0] - _caps[1]);
  Real radius = _radius;

  glPushMatrix();
    VEC3F axis;
    Real angle;
    MATRIX_UTIL::quaternionToAxisAngle(_cylinderRotation, axis, angle);
    // _cylinderRotation.axisAngle(axis, angle);
    glTranslatef(_cylinderTranslation[0], _cylinderTranslation[1], _cylinderTranslation[2]);
    glRotatef(angle, axis[0], axis[1], axis[2]);

    Real cylHalfHeight = length / 2.0;

    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i < 33; i++)
    {
      Real angle = i / 32.0 * 2.0 * M_PI;
      Real ca = cos(angle);
      Real sa = sin(angle);

      VEC3F normal(ca, sa, 0);
      normal = _cylinderRotation.toRotationMatrix() * normal;
      glNormal3f(normal[0], normal[1], normal[2]);
 
      glVertex3f(radius * ca, cylHalfHeight, radius * sa);
      glVertex3f(radius * ca, -cylHalfHeight, radius * sa);
    }
    glEnd();

    // draw a line towards one of the end caps
    glDisable(GL_DEPTH_TEST);
    glColor4f(0,10,0,1);
    glLineWidth(10);
    glBegin(GL_LINES);
      glVertex3f(0, 0, 0);
      glVertex3f(0, cylHalfHeight, 0);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  glPopMatrix();
}

//////////////////////////////////////////////////////////////////////
// return the bounding box dimensions
//////////////////////////////////////////////////////////////////////
void CYLINDER::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  // get the eight corner of the original bounding box
  VEC3F corners[8];

  Real distance = sqrt(_radius * _radius + _radius * _radius);

  Real length = fabs(_caps[0] - _caps[1]);

  corners[0] = VEC3F(distance, length * 0.5, distance); 
  corners[1] = VEC3F(-distance, length * 0.5, distance); 
  corners[2] = VEC3F(distance, length * 0.5, -distance); 
  corners[3] = VEC3F(-distance, length * 0.5, -distance);

  corners[4] = VEC3F(distance, -length * 0.5, distance); 
  corners[5] = VEC3F(-distance, -length * 0.5, distance); 
  corners[6] = VEC3F(distance, -length * 0.5, -distance); 
  corners[7] = VEC3F(-distance, -length * 0.5, -distance);

  // transform the corners
  MATRIX3 rotation = _cylinderRotation.toRotationMatrix();
  for (int x = 0; x < 8; x++)
    corners[x] =  rotation * corners[x] + _cylinderTranslation;

  // find the bounds
  mins = corners[0];
  maxs = corners[0];

  for (int x = 1; x < 8; x++)
    for (int y = 0; y < 3; y++)
    {
      if (corners[x][y] < mins[y])
        mins[y] = corners[x][y];
      if (corners[x][y] > maxs[y])
        maxs[y] = corners[x][y];
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CYLINDER::setLength(Real length)
{
  Real middle = (_caps[1] + _caps[0]) * 0.5;

  _caps[0] = middle - length * 0.5;
  _caps[1] = middle + length * 0.5;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CYLINDER::endPoints(VEC3F& leftCap, VEC3F& rightCap)
{
  VEC3F& translation = _cylinderTranslation;
  MATRIX3 rotation = _cylinderRotation.toRotationMatrix();
  Real length = getLength();
  VEC3F beginVertex(0,0, length * 0.5 + _radius);
  VEC3F endVertex(0,0, -length * 0.5 - _radius);

  // apply the ODE-style rotation first
  rotation = rotation * Eigen::AngleAxisd(0.5 * M_PI, VEC3F::UnitX());

  leftCap  = rotation * beginVertex + translation;
  rightCap = rotation * endVertex + translation;
}
