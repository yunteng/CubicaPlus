// SETTINGS.h: Project-wide options set in one place
//
//////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

#define USING_OSX __APPLE__

#define USING_GLVU 1

// select single or double precision
#define DOUBLE_PRECISION 1

#define USING_OPENMP 1
#define USING_SUBSPACE_OPENMP 1

typedef Eigen::Vector3i VEC3I;
#if DOUBLE_PRECISION
typedef Eigen::Vector3d VEC3F;
typedef Eigen::Vector2d VEC2F;
typedef Eigen::VectorXd VECTOR;
typedef Eigen::ArrayXd ARRAY;
typedef Eigen::MatrixXd MATRIX;
typedef double Real;
typedef Eigen::Matrix3d MATRIX3;
typedef Eigen::Matrix4d MATRIX4;
typedef Eigen::Matrix<double, 3, 4, Eigen::ColMajor> MATRIX3x4;
typedef Eigen::Matrix<double, 9, 9, Eigen::RowMajor> MATRIX9;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Quaternion<double> QUATERNION;
typedef Eigen::Triplet<double> TRIPLET;
// typedef Eigen::DiagonalMatrix<double> DiagMat;
#else
typedef Eigen::Vector2f VEC2F;
typedef Eigen::Vector3f VEC3F;
typedef Eigen::VectorXf VECTOR;
typedef Eigen::ArrayXf ARRAY;
typedef Eigen::MatrixXf MATRIX;
typedef float Real;
typedef Eigen::Matrix3f MATRIX3;
typedef Eigen::Matrix4f MATRIX4;
typedef Eigen::Matrix<float, 9, 9, Eigen::RowMajor> MATRIX9;
typedef Eigen::Quaternion<float> QUATERNION;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<float> TRIPLET;
// typedef Eigen::DiagonalMatrix<float> DiagMat;
#endif

// Any squared distance measured between two points will be at least this much.
const double FORCE_GENERAL_POSITION_EPSILON_SQR = 1e-8;

// turn off asserts?
#define NDEBUG

#endif
