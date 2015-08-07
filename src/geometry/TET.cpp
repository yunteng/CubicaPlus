/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <geometry/TET.h>
#include <util/MATRIX_UTIL.h>
TET::TET(VEC3F& v0, VEC3F& v1, VEC3F& v2, VEC3F& v3):
  _materialIndex(0)
{
  vertices[0] = &v0;
  vertices[1] = &v1;
  vertices[2] = &v2;
  vertices[3] = &v3;
  for(int x = 0; x < 4; x++)
    _restVertices[x] = NULL;
}

TET::TET(const TET& tet)
{
  vertices[0] = tet.vertices[0];
  vertices[1] = tet.vertices[1];
  vertices[2] = tet.vertices[2];
  vertices[3] = tet.vertices[3];
  for(int x = 0; x < 4; x++)
    _restVertices[x] = tet._restVertices[x];
  _materialIndex = tet._materialIndex;
}

TET::~TET()
{
  
}

void TET::init()
{
  if(_restVertices[0] == NULL) return;

  vector<TRIANGLE> faces;
  for(int x = 0; x < 4; x++)
    faces.push_back(face(x));

  VEC3F row[3];
  row[0] = *(_restVertices[1]) - *(_restVertices[0]);
  row[1] = *(_restVertices[2]) - *(_restVertices[0]),
  row[2] = *(_restVertices[3]) - *(_restVertices[0]);

  _DmInv.setZero();

  _DmInv << row[0][0], row[1][0], row[2][0], 
            row[0][1], row[1][1], row[2][1], 
            row[0][2], row[1][2], row[2][2];

  _DmInv = _DmInv.inverse().eval();

  // calculate the area vectors
  // v0 is incident on faces (0,1,2)
  _b[0] = faces[0].normal() * faces[0].area() +
          faces[1].normal() * faces[1].area() +
          faces[2].normal() * faces[2].area();
  
  // v1 is incident on faces (0,1,3)
  _b[1] = faces[0].normal() * faces[0].area() +
          faces[1].normal() * faces[1].area() +
          faces[3].normal() * faces[3].area();
  
  // v2 is incident on faces (1,2,3)
  _b[2] = faces[1].normal() * faces[1].area() +
          faces[2].normal() * faces[2].area() +
          faces[3].normal() * faces[3].area();
  
  // v3 is incident on faces (0,2,3)
  _b[3] = faces[0].normal() * faces[0].area() +
          faces[2].normal() * faces[2].area() +
          faces[3].normal() * faces[3].area();
  _b[0] *= -1.0f / 3.0f;
  _b[1] *= -1.0f / 3.0f;
  _b[2] *= -1.0f / 3.0f;
  _b[3] *= -1.0f / 3.0f;

  for(int y = 0; y < 4; y++)
    _bMat.col(y) = _b[y];


  _restVolume = fabs(row[0].dot(row[1].cross(row[2])) / 6.0);

  computePFPu();
}

//////////////////////////////////////////////////////////////////////
// compute the tet volume
//////////////////////////////////////////////////////////////////////
Real TET::volume()
{
  // formula for a tet volume with vertices (a,b,c,d) is:
  // |(a - d) dot ((b - d) cross (c - d))| / 6
  VEC3F a = (*vertices[1]) - (*vertices[0]);
  VEC3F b = (*vertices[2]) - (*vertices[0]);
  VEC3F c = (*vertices[3]) - (*vertices[0]);

  return fabs(a.dot(b.cross(c)) / 6.0);
}

VEC3F TET::center()
{
  VEC3F sum;
  sum += *(vertices[0]);
  sum += *(vertices[1]);
  sum += *(vertices[2]);
  sum += *(vertices[3]);
  sum *= 0.25;
  return sum;
}

VEC3F TET::restCenter()
{
  VEC3F sum;
  sum += *(_restVertices[0]);
  sum += *(_restVertices[1]);
  sum += *(_restVertices[2]);
  sum += *(_restVertices[3]);
  sum *= 0.25;
  return sum;
}

TRIANGLE TET::face(int x)
{
  if (x == 0) return TRIANGLE(*vertices[0], *vertices[1], *vertices[3]);
  if (x == 1) return TRIANGLE(*vertices[0], *vertices[2], *vertices[1]);
  if (x == 2) return TRIANGLE(*vertices[3], *vertices[2], *vertices[0]);
  return TRIANGLE(*vertices[1], *vertices[2], *vertices[3]);
}

MATRIX3 TET::F() {
  // Following Teran's notation - spatial tensor
  VEC3F row[3];
  row[0] = *(vertices[1]) - *(vertices[0]);
  row[1] = *(vertices[2]) - *(vertices[0]),
  row[2] = *(vertices[3]) - *(vertices[0]);

  MATRIX3 Ds;
  Ds << row[0][0], row[1][0], row[2][0], 
        row[0][1], row[1][1], row[2][1], 
        row[0][2], row[1][2], row[2][2];

  return Ds * _DmInv;
}
MATRIX3 TET::G()
{
  VEC3F row[3];
  row[0] = *(vertices[0]) - *(vertices[3]);
  row[1] = *(vertices[1]) - *(vertices[3]),
  row[2] = *(vertices[2]) - *(vertices[3]);

  MATRIX3 T;
  T << row[0][0], row[1][0], row[2][0], 
       row[0][1], row[1][1], row[2][1], 
       row[0][2], row[1][2], row[2][2];
  return T.transpose().inverse();
}
//////////////////////////////////////////////////////////////////////
// calculate the barycentric coordinates of a vertex, if rest = true, 
//this is done in rest position, otherwise in current position
//////////////////////////////////////////////////////////////////////
bool TET::baryCenter(const VEC3F& vert, QUATERNION& lambda){
  VEC3F row[3];
  row[0] = *(vertices[0]) - *(vertices[3]);
  row[1] = *(vertices[1]) - *(vertices[3]),
  row[2] = *(vertices[2]) - *(vertices[3]);

  MATRIX3 T;
  T << row[0][0], row[1][0], row[2][0], 
       row[0][1], row[1][1], row[2][1], 
       row[0][2], row[1][2], row[2][2];


  VEC3F lambda0_3 = T.inverse() * (vert - *(vertices[3]));
  lambda.x() = lambda0_3[0];
  lambda.y() = lambda0_3[1];
  lambda.z() = lambda0_3[2];
  lambda.w() = 1.0 - lambda.x() - lambda.y() - lambda.z();

  if(!(
      0 <= lambda.x() && lambda.x() <= 1
   && 0 <= lambda.y() && lambda.y() <= 1
   && 0 <= lambda.z() && lambda.z() <= 1
   && 0 <= lambda.w() && lambda.w() <= 1) )
    return false;
  return true;
  
}

bool TET::restBaryCenter(const VEC3F& vert, QUATERNION& lambda)
{
  VEC3F row[3];
  row[0] = *(_restVertices[0]) - *(_restVertices[3]);
  row[1] = *(_restVertices[1]) - *(_restVertices[3]),
  row[2] = *(_restVertices[2]) - *(_restVertices[3]);

  MATRIX3 T;
  T << row[0][0], row[1][0], row[2][0], 
       row[0][1], row[1][1], row[2][1], 
       row[0][2], row[1][2], row[2][2];

  VEC3F lambda0_3 = T.inverse() * (vert - *(_restVertices[3]));
  lambda.x() = lambda0_3[0];
  lambda.y() = lambda0_3[1];
  lambda.z() = lambda0_3[2];
  lambda.w() = 1.0 - lambda.x() - lambda.y() - lambda.z();

  if(!(
      0 <= lambda.x() && lambda.x() <= 1
   && 0 <= lambda.y() && lambda.y() <= 1
   && 0 <= lambda.z() && lambda.z() <= 1
   && 0 <= lambda.w() && lambda.w() <= 1) )
    return false;
  return true;
}

void TET::drawFaces()
{
  for(int x = 0; x < 4; x++)
    face(x).draw();
}

//////////////////////////////////////////////////////////////////////
// This computes \partial{F} / \partial{u} so that if a derivative
// is taken with respect to F (the deformation gradient), it can
// be converted via the chain rule to a derivative wrt u. The symbolic
// derivative wrt F is less memory intensive since it comprises a
// smaller number of symbols overall.
//
// autogenerated from:
//
/*
function code = codegen_DFDu()
    a = sym('a'); b = sym('b'); c = sym('c'); d = sym('d');
    e = sym('e'); f = sym('f'); g = sym('g'); h = sym('h');
    i = sym('i'); j = sym('j'); k = sym('k'); l = sym('l');    
    m = sym('m'); n = sym('n'); o = sym('o'); p = sym('p');
    q = sym('q'); r = sym('r'); s = sym('s'); t = sym('t');
    v = sym('u');    
    n0 = [a b c];
    n1 = [d e f];
    n2 = [g h i];
    n3 = [j k l];
    u = mytrans([a b c d e f g h i j k l]);
    x = mytrans([
        n1-n0
        n2-n0
        n3-n0
        ]);
    Xinv = [
        m n o
        p q r
        s t v
        ];

    F = x * Xinv;
    F = reshape(F, 9, 1);
    
    pF_pu = sym('a') * ones(9, 12);
    for i = 1:9
        for j = 1:12
            pF_pu(i,j) = diff( F(i), u(j) );
        end
    end    
    code = maple('codegen[C]', pF_pu, 'optimized');
end
*/
void TET::computePFPu()
{
  const MATRIX3& matInv = _DmInv;
  const double m = matInv(0,0);
  const double n = matInv(0,1);
  const double o = matInv(0,2);
  const double p = matInv(1,0);
  const double q = matInv(1,1);
  const double r = matInv(1,2);
  const double s = matInv(2,0);
  const double t = matInv(2,1);
  const double u = matInv(2,2);

  const double t1 = -m-p-s;
  const double t2 = -n-q-t;
  const double t3 = -o-r-u;

  _PFPu.resize(9, 12);
  _PFPu.setZero();
  _PFPu(0,0) = t1;
  _PFPu(0,1) = 0.0;
  _PFPu(0,2) = 0.0;
  _PFPu(0,3) = m;
  _PFPu(0,4) = 0.0;
  _PFPu(0,5) = 0.0;
  _PFPu(0,6) = p;
  _PFPu(0,7) = 0.0;
  _PFPu(0,8) = 0.0;
  _PFPu(0,9) = s;
  _PFPu(0,10) = 0.0;
  _PFPu(0,11) = 0.0;
  _PFPu(1,0) = 0.0;
  _PFPu(1,1) = t1;
  _PFPu(1,2) = 0.0;
  _PFPu(1,3) = 0.0;
  _PFPu(1,4) = m;
  _PFPu(1,5) = 0.0;
  _PFPu(1,6) = 0.0;
  _PFPu(1,7) = p;
  _PFPu(1,8) = 0.0;
  _PFPu(1,9) = 0.0;
  _PFPu(1,10) = s;
  _PFPu(1,11) = 0.0;
  _PFPu(2,0) = 0.0;
  _PFPu(2,1) = 0.0;
  _PFPu(2,2) = t1;
  _PFPu(2,3) = 0.0;
  _PFPu(2,4) = 0.0;
  _PFPu(2,5) = m;
  _PFPu(2,6) = 0.0;
  _PFPu(2,7) = 0.0;
  _PFPu(2,8) = p;
  _PFPu(2,9) = 0.0;
  _PFPu(2,10) = 0.0;
  _PFPu(2,11) = s;
  _PFPu(3,0) = t2;
  _PFPu(3,1) = 0.0;
  _PFPu(3,2) = 0.0;
  _PFPu(3,3) = n;
  _PFPu(3,4) = 0.0;
  _PFPu(3,5) = 0.0;
  _PFPu(3,6) = q;
  _PFPu(3,7) = 0.0;
  _PFPu(3,8) = 0.0;
  _PFPu(3,9) = t;
  _PFPu(3,10) = 0.0;
  _PFPu(3,11) = 0.0;
  _PFPu(4,0) = 0.0;
  _PFPu(4,1) = t2;
  _PFPu(4,2) = 0.0;
  _PFPu(4,3) = 0.0;
  _PFPu(4,4) = n;
  _PFPu(4,5) = 0.0;
  _PFPu(4,6) = 0.0;
  _PFPu(4,7) = q;
  _PFPu(4,8) = 0.0;
  _PFPu(4,9) = 0.0;
  _PFPu(4,10) = t;
  _PFPu(4,11) = 0.0;
  _PFPu(5,0) = 0.0;
  _PFPu(5,1) = 0.0;
  _PFPu(5,2) = t2;
  _PFPu(5,3) = 0.0;
  _PFPu(5,4) = 0.0;
  _PFPu(5,5) = n;
  _PFPu(5,6) = 0.0;
  _PFPu(5,7) = 0.0;
  _PFPu(5,8) = q;
  _PFPu(5,9) = 0.0;
  _PFPu(5,10) = 0.0;
  _PFPu(5,11) = t;
  _PFPu(6,0) = t3;
  _PFPu(6,1) = 0.0;
  _PFPu(6,2) = 0.0;
  _PFPu(6,3) = o;
  _PFPu(6,4) = 0.0;
  _PFPu(6,5) = 0.0;
  _PFPu(6,6) = r;
  _PFPu(6,7) = 0.0;
  _PFPu(6,8) = 0.0;
  _PFPu(6,9) = u;
  _PFPu(6,10) = 0.0;
  _PFPu(6,11) = 0.0;
  _PFPu(7,0) = 0.0;
  _PFPu(7,1) = t3;
  _PFPu(7,2) = 0.0;
  _PFPu(7,3) = 0.0;
  _PFPu(7,4) = o;
  _PFPu(7,5) = 0.0;
  _PFPu(7,6) = 0.0;
  _PFPu(7,7) = r;
  _PFPu(7,8) = 0.0;
  _PFPu(7,9) = 0.0;
  _PFPu(7,10) = u;
  _PFPu(7,11) = 0.0;
  _PFPu(8,0) = 0.0;
  _PFPu(8,1) = 0.0;
  _PFPu(8,2) = t3;
  _PFPu(8,3) = 0.0;
  _PFPu(8,4) = 0.0;
  _PFPu(8,5) = o;
  _PFPu(8,6) = 0.0;
  _PFPu(8,7) = 0.0;
  _PFPu(8,8) = r;
  _PFPu(8,9) = 0.0;
  _PFPu(8,10) = 0.0;
  _PFPu(8,11) = u;
}
