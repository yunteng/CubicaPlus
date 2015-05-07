#include <geometry/TRIANGLE.h>
#include <iostream>
#include <algorithm>

TRIANGLE::TRIANGLE(VEC3F& v0, VEC3F& v1, VEC3F& v2, VEC3F color):
  _color(color)
{
  vertices[0] = &v0;
  vertices[1] = &v1;
  vertices[2] = &v2;
  _normal = (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]);
  _area = 0.5 * _normal.norm();
  _normal.normalize();
}

TRIANGLE::TRIANGLE(const TRIANGLE& triangle)
{
  _color = triangle._color;
  vertices[0] = triangle.vertices[0];
  vertices[1] = triangle.vertices[1];
  vertices[2] = triangle.vertices[2];
  _normal = triangle._normal;
  _area = triangle._area;
}

void TRIANGLE::draw()
{
  _normal.normalize();
  glBegin(GL_TRIANGLES);
    glNormal3f(_normal[0], _normal[1], _normal[2]);
    VEC3F& v0 = *vertices[0];
    VEC3F& v1 = *vertices[1];
    VEC3F& v2 = *vertices[2];

    glVertex3f(v0[0], v0[1], v0[2]);
    glVertex3f(v1[0], v1[1], v1[2]);
    glVertex3f(v2[0], v2[1], v2[2]);
  glEnd();
}

void TRIANGLE::draw(const vector<VEC3F>& colors)
{
  glBegin(GL_TRIANGLES);
    glNormal3f(_normal[0], _normal[1], _normal[2]);
    VEC3F& v0 = *vertices[0];
    VEC3F& v1 = *vertices[1];
    VEC3F& v2 = *vertices[2];

    glColor3f(colors[0][0], colors[0][1], colors[0][2]);
    glVertex3f(v0[0], v0[1], v0[2]);
    glColor3f(colors[1][0], colors[1][1], colors[1][2]);
    glVertex3f(v1[0], v1[1], v1[2]);
    glColor3f(colors[2][0], colors[2][1], colors[2][2]);
    glVertex3f(v2[0], v2[1], v2[2]);
  glEnd();
}
//////////////////////////////////////////////////////////////////////
// for equality, sort the indices first and then determine
// if they are all the same
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::operator==(const TRIANGLE& RHS) const {
  vector<VEC3F*> copyLHS(vertices, vertices + 3);;
  vector<VEC3F*> copyRHS(RHS.vertices, RHS.vertices + 3);
  sort(copyLHS.begin(), copyLHS.end());
  sort(copyRHS.begin(), copyRHS.end());

  for (unsigned int x =0; x < copyLHS.size(); x++)
    if (copyLHS[x] != copyRHS[x])
      return false;

  return true;
}
VEC3F TRIANGLE::center()
{
  VEC3F sum;
  sum.setZero();
  sum += *(vertices[0]);
  sum += *(vertices[1]);
  sum += *(vertices[2]);
  sum *= 1.0 / 3.0;
  return sum;
}
void TRIANGLE::boundingBox(VEC3F& mins, VEC3F& maxs)
{
  mins = *vertices[0];
  maxs = *vertices[0];

  for (int x = 1; x < 3; x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (*vertices[x])[y] < mins[y] ? (*vertices[x])[y] : mins[y];
      maxs[y] = (*vertices[x])[y] > maxs[y] ? (*vertices[x])[y] : maxs[y];
    }
}
bool TRIANGLE::baryCenter(const VEC3F& point, VEC3F& lambda)
{
  VEC3F& a = *vertices[0];
  VEC3F& b = *vertices[1];
  VEC3F& c = *vertices[2];

  VEC3F vecA = (b-a).cross(c - a);
  VEC3F vecAa = (b - point).cross(c - point);
  VEC3F vecAb = (point - a).cross(c - a);
  VEC3F vecAc = (b - a).cross(point - a);

  if(vecAa.dot(vecA) < 0 || vecAb.dot(vecA) < 0 || vecAc.dot(vecA) < 0){
    return false;
  }
  
  Real A = vecA.norm();
  Real Aa = vecAa.norm();
  Real Ab = vecAb.norm();
  Real Ac = vecAc.norm();
  lambda[0] = Aa / A;
  lambda[1] = Ab / A;
  lambda[2] = Ac / A;
  Real length = lambda[0] + lambda[1] + lambda[2];
  if(abs(length - 1) <= 0.1){
    lambda /= length;
    return true;
  }
  return false;
}
VEC3F TRIANGLE::interpPosition(VEC3F& lambda){
  return lambda[0] * (*vertices[0]) + lambda[1] * (*vertices[1]) + lambda[2] * (*vertices[2]);
}
//////////////////////////////////////////////////////////////////////
// intersect with line segment
//////////////////////////////////////////////////////////////////////
bool TRIANGLE::intersects(const VEC3F& start, const VEC3F& end)
{
  VEC3F& a = *vertices[0];
  VEC3F& b = *vertices[1];
  VEC3F& c = *vertices[2];

  VEC3F geometricNormal = (b - a).cross(c - a);
  geometricNormal.normalize();

  VEC3F direction = end - start;
  Real length = direction.norm();
  direction.normalize();
  VEC3F diff = a - start;
  double denom = direction.dot(geometricNormal);

  // catch divide by zero
  if (fabs(denom) <= 0.0) return false;

  double t = diff.dot(geometricNormal) / denom;
  if (t < 0) return false;

  VEC3F h = start + (direction * t);

  VEC3F test = (b - a).cross(h - a);
  if (test.dot(geometricNormal) < 0) return false; 
  test = (c - b).cross(h - b);
  if (test.dot(geometricNormal) < 0) return false; 
  test = (a - c).cross(h - c);
  if (test.dot(geometricNormal) < 0) return false; 

  return (t <= length);
}
VEC3F TRIANGLE::projection(const VEC3F& point)
{
  VEC3F v = point - *vertices[0];
  _normal = normal();
  return point - v.dot(_normal) * _normal;
}
