/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#include <SETTINGS.h>
#include <Eigen/SVD>
#include <linearalgebra/COO_MATRIX.h>

class MATRIX_UTIL
{
public:
  static MATRIX3 repackVec9(const VECTOR& F)
  {
    MATRIX3 final;
    final(0,0) = F(0);
    final(1,0) = F(1);
    final(2,0) = F(2);
    final(0,1) = F(3);
    final(1,1) = F(4);
    final(2,1) = F(5);
    final(0,2) = F(6);
    final(1,2) = F(7);
    final(2,2) = F(8);

    return final;
  };

  static MATRIX3 repackVec9(const Real* F)
  {
    MATRIX3 final;
    final(0,0) = F[0];
    final(1,0) = F[1];
    final(2,0) = F[2];
    final(0,1) = F[3];
    final(1,1) = F[4];
    final(2,1) = F[5];
    final(0,2) = F[6];
    final(1,2) = F[7];
    final(2,2) = F[8];

    return final;
  };

  static void repackVec9(const VECTOR& F, MATRIX3& final)
  {
    final(0,0) = F(0);
    final(1,0) = F(1);
    final(2,0) = F(2);
    final(0,1) = F(3);
    final(1,1) = F(4);
    final(2,1) = F(5);
    final(0,2) = F(6);
    final(1,2) = F(7);
    final(2,2) = F(8);
  };

  static VECTOR flatternMatrix3(const MATRIX3 F)
  {
    VECTOR final(9);
    final(0) = F(0,0);
    final(1) = F(1,0);
    final(2) = F(2,0);
    final(3) = F(0,1);
    final(4) = F(1,1);
    final(5) = F(2,1);
    final(6) = F(0,2);
    final(7) = F(1,2);
    final(8) = F(2,2);

    return final;
  };

  static void flatternMatrix3(const MATRIX3& F, VECTOR& output)
  {
    if(output.size() != 9 )
      output.resize(9);

    output(0) = F(0,0);
    output(1) = F(1,0);
    output(2) = F(2,0);
    output(3) = F(0,1);
    output(4) = F(1,1);
    output(5) = F(2,1);
    output(6) = F(0,2);
    output(7) = F(1,2);
    output(8) = F(2,2);
  };

  static MATRIX3 cross(const VEC3F& vec)
  {
    MATRIX3 final;
    final(0,0) = 0.0;
    final(1,0) = vec[2];
    final(2,0) = -vec[1];

    final(0,1) = -vec[2];
    final(1,1) = 0.0;
    final(2,1) = vec[0];

    final(0,2) = vec[1];
    final(1,2) = -vec[0];
    final(2,2) = 0.0;

    return final;
  }
  static MATRIX3 outer_product(const VEC3F& v)
  {
    MATRIX3 A;
    const Real& x=v[0], y=v[1], z=v[2];

    A(0,0) = x*x;  A(0,1) = x*y;  A(0,2) = x*z;
    A(1,0)=A(0,1); A(1,1) = y*y;  A(1,2) = y*z;
    A(2,0)=A(0,2); A(2,1)=A(1,2); A(2,2) = z*z;

    return A;
  }
  static void vectorToMatrix(const vector<VECTOR>& vectors, MATRIX& mat, bool asColumns)
  {
    // each vector is a column of the matrix
    if(asColumns)
      mat.resize(vectors[0].size(), vectors.size());
    else
      mat.resize(vectors.size(), vectors[0].size());

    for(unsigned int x = 0; x < vectors.size(); x++){
      if(asColumns){
        mat.col(x) = vectors[x];
      }else{
        mat.row(x) = vectors[x];
      }
    }
  }
  static void matrixToVectors(const MATRIX& mat, vector<VECTOR>& vectors, bool asColumns)
  {
    if(asColumns){
      vectors.resize(mat.cols());
      for(int x = 0; x < mat.cols(); x++)
        vectors[x] = mat.col(x);
    }else{
      vectors.resize(mat.rows());
      for(int x = 0; x < mat.rows(); x++)
        vectors[x] = mat.row(x);
    }

  }
  static void pcaEigen(const vector<VECTOR>& data, bool shift, MATRIX& components, VECTOR& values)
  {
    MATRIX mat;
    vectorToMatrix(data, mat, true);
    pcaEigen(mat, shift, components, values);
  }
  static void pcaEigen(MATRIX& data, bool shift, MATRIX& components, VECTOR& values)
  {
    int rows = data.rows();
    int cols = data.cols();
    if(shift){
      for(int x = 0; x < rows; x++){
        Real mean = data.row(x).sum() / cols;
        data.row(x) = (data.row(x).array() - mean).matrix();
      }
      data *= 1.0 / sqrt(cols - 1);
    }
    Eigen::JacobiSVD<MATRIX> eigenSystem(data, Eigen::ComputeThinU);
    const VECTOR& singularValues = eigenSystem.singularValues();

    components = eigenSystem.matrixU();

    values.resize(components.cols());
    for(unsigned int x = 0; x < components.cols(); x++)
      values[x] = singularValues[x] * singularValues[x];
  }

  static void verifyTruncatedPCA(const MATRIX& originalData, const MATRIX& U)
  {
    Real meanDiffNorm = 0;
    Real maxDiffNorm = 0;
    Real meanRelativeDiffNorm = 0;
    Real maxRelativeDiffNorm = 0;
    vector<Real> diffNorms;
    vector<Real> relativeDiffNorms;

    Real maxRelativeActualNorm;
    for (int i = 0; i < originalData.cols(); i++)
    {
      VECTOR testColumn = originalData.col(i);

      VECTOR testProjection = U * (U.transpose() * testColumn);

      VECTOR diff = testColumn - testProjection;
      
      diffNorms.push_back(diff.norm());

      meanDiffNorm += diff.norm();

      Real relativeError = diff.norm() / testColumn.norm();
      meanRelativeDiffNorm += relativeError;

      relativeDiffNorms.push_back(relativeError);

      if (diff.norm() > maxDiffNorm)
        maxDiffNorm = diff.norm();

      if (relativeError > maxRelativeDiffNorm)
      {
        maxRelativeDiffNorm = relativeError;
        maxRelativeActualNorm = maxDiffNorm;
      }
    }

    cout << "======================================================" << endl;
    cout << " PCA test results " << endl;
    cout << "  Mean absolute projection error: " << meanDiffNorm / originalData.cols() << endl;
    cout << "  Max absolute projection error: " << maxDiffNorm << endl ;
    cout << "  Mean relative projection error: " << meanRelativeDiffNorm / originalData.cols() << endl;
    cout << "  Max relative projection error: " << maxRelativeDiffNorm << endl;
    cout << "  Actual norm of max relative error sample: " << maxRelativeActualNorm << endl;

    cout << "======================================================" << endl;
  }

  static void orthogonalize(MATRIX& mat, int c = 0)
  {
    if(c >= mat.cols())
      return;

    for (int x = c; x < mat.cols(); x++)
    {
      // VECTOR& currentColumn = mat.col(x);

      // project out other components
      for (int y = 0; y < x; y++)
      {
        // VECTOR& previousColumn = mat.col(y);
        Real dot = mat.col(y).dot(mat.col(x));
        mat.col(x) -= dot * mat.col(y);
      }

      Real norm2 = mat.col(x).norm();
      // assert(norm2 != 0);
      if(norm2 < 1e-6){
        cout << "zero column " << x << endl;
      }else
        mat.col(x) *= 1.0 / norm2;
    }
  }

  static void eigensystem3x3(const MATRIX& mat, VEC3F& eigenvalues, MATRIX3& eigenvectors)
  {
    register int i, j, k, l;

    //float TOL = 1e-8f;
    //int MAX_SWEEPS = 500;
    float TOL = 1e-3f;
    int MAX_SWEEPS = 50;
    unsigned int n = 3;

    float a[9];
    float d[3];
    float v[9];
    i = 0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++, i++)
      {
        a[i] = mat(x,y);
        v[i] = (x == y) ? 1.0 : 0.0;
      }
    
    
    float onorm, dnorm;
    float b, dma, q, t, c, s;
    float atemp, vtemp, dtemp;

    // Set v to the identity matrix, set the vector d to contain the
    // diagonal elements of the matrix a
    d[0] = a[0];
    d[1] = a[4];
    d[2] = a[8];

    for (l = 1; l <= MAX_SWEEPS; l++)
    {
      // Set dnorm to be the maximum norm of the diagonal elements, set
      // onorm to the maximum norm of the off-diagonal elements
      
      dnorm = (float)fabs(d[0]) + (float)fabs(d[1]) + (float)fabs(d[2]);
      onorm = (float)fabs(a[1]) + (float)fabs(a[2]) + (float)fabs(a[5]);
      // Normal end point of this algorithm.
      if((onorm/dnorm) <= TOL)
        goto Exit_now;

      for (j = 1; j < static_cast<int>(n); j++)
      {
        for (i = 0; i <= j - 1; i++)
        {

          b = a[n*i+j];
          if(fabs(b) > 0.0f)
          {
            dma = d[j] - d[i];
            if((fabs(dma) + fabs(b)) <= fabs(dma))
              t = b / dma;
            else
            {
              q = 0.5f * dma / b;
              t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
              if (q < 0.0)
                t = -t;
            }

            c = 1.0f/(float)sqrt(t*t + 1.0f);
            s = t * c;
            a[n*i+j] = 0.0f;

            for (k = 0; k <= i-1; k++)
            {
              atemp = c * a[n*k+i] - s * a[n*k+j];
              a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
              a[n*k+i] = atemp;
            }

            for (k = i+1; k <= j-1; k++)
            {
              atemp = c * a[n*i+k] - s * a[n*k+j];
              a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
              a[n*i+k] = atemp;
            }

            for (k = j+1; k < static_cast<int>(n); k++)
            {
              atemp = c * a[n*i+k] - s * a[n*j+k];
              a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
              a[n*i+k] = atemp;
            }

            for (k = 0; k < static_cast<int>(n); k++)
            {
              vtemp = c * v[n*k+i] - s * v[n*k+j];
              v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
              v[n*k+i] = vtemp;
            }

            dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
            d[j] = s*s*d[i] + c*c*d[j] + 2.0f*c*s*b;
            d[i] = dtemp;
          }
        }
      }
    }

  Exit_now:
    for (int x = 0; x < 3; x++)
      eigenvalues[x] = d[x];

    i = 0;
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++, i++)
        eigenvectors(x,y) = v[i];

    return;
  }

  static void quaternionToAxisAngle(const QUATERNION& quat, VEC3F& axis, Real& angle)
  {
    Real _w = quat.w();
    Real _x = quat.x();
    Real _y = quat.y();
    Real _z = quat.z();
    if ((_w >= ((Real)1)) || (_w <= (Real)(-1)))
    {
      // identity; this check is necessary to avoid problems with acos if s is 1 + eps
      angle = 0;
      axis[0] = 1;
      axis[1] = 0;
      axis[2] = 0;
      //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      //cout << "Quaternion is not a rotation! " << endl;
      return;
    }

    angle = 2.0 * acos(_w);
    Real sin2 = _x*_x + _y*_y + _z*_z; //sin^2(*angle / 2.0)

    if (sin2 == 0)
    {
      // identity rotation; angle is zero, any axis is equally good
      axis[0] = 1;
      axis[1] = 0;
      axis[2] = 0;
    }
    else
    {
      Real inv = 1.0 / sqrt(sin2); // note: *angle / 2.0 is on [0,pi], so sin(*angle / 2.0) >= 0, and therefore the sign of sqrt can be safely taken positive
      axis[0] = _x * inv;
      axis[1] = _y * inv;
      axis[2] = _z * inv;
    }

    angle = angle / (2.0 * M_PI) * 360.0f;
  }

  // static void SpSymMatToUmSymMat(const SpMat& input, ARumSymMatrix<Real>& output)
  // {
  //  int n = input.rows();
  //  int nnz = (input.nonZeros() + n) / 2;
  //  Real 
  //  ARumSymMatrix(int np, int nnzp, ARTYPE* ap, int* irowp,
 //                int* pcolp, char uplop = 'L', double thresholdp = 0.1,
 //                int fillinp = 9, bool reducible = true, bool check = true);
  // }
};
inline QUATERNION operator-(const QUATERNION& a, const QUATERNION& b)
{
  return QUATERNION(a.w() - b.w(),
                    a.x() - b.x(),
                    a.y() - b.y(),
                    a.z() - b.z());
}
inline QUATERNION operator*(const QUATERNION& q, const Real s)
{
  return QUATERNION(q.w() * s, q.x() * s, q.y() * s, q.z() * s);
}
inline QUATERNION operator+(const QUATERNION& a, const QUATERNION& b)
{
  return QUATERNION(a.w() + b.w(), 
                    a.x() + b.x(),
                    a.y() + b.y(),
                    a.z() + b.z());
}

inline std::istream &operator>>(std::istream &in, VEC3F& v)
{ return in >> v[0] >> v[1] >> v[2]; }

inline std::istream &operator>>(std::istream &in, VEC2F& v)
{ return in >> v[0] >> v[1]; }

inline VEC3F cross(const VEC3F& a, const VEC3F& b){
  return a.cross(b);
}
inline void unitize(VEC3F a)
{
  a.normalize();
}
inline Real norm(const VEC3F& a)
{
  return a.norm();
}
inline Real norm2(const VEC3F& a)
{
  return a.dot(a);
}

inline VEC3F operator^(const VEC3F& a, const VEC3F& b)
{
  return a.cross(b);
}
inline MATRIX3 operator^(const MATRIX3& n, const MATRIX3& m)
{
  return n.transpose() * m;
}
inline Real trace(const MATRIX3& mat)
{
  return mat.trace();
}
inline int lowerNNz(const MATRIX& mat){
  return mat.rows() * (mat.rows() + 1) * 0.5;
}
inline int lowerNNz(const COO_MATRIX& mat){
  return mat.lowerNNz();
}
// inline MATRIX3& operator+=(MATRIX3& m, Real s)
// { m(0, 0) += s; m(0, 1) += s; m(0, 2) += s; 
//  m(1, 0) += s; m(1, 1) += s; m(1, 2) += s;
//  m(2, 0) += s; m(2, 1) += s; m(2, 2) += s;
//  return m; 
// }

#endif
