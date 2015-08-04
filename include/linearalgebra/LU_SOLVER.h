/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef LU_SOLVER_H
#define LU_SOLVER_H

#include <SETTINGS.h>

#if USING_OSX
#include <Accelerate/Accelerate.h>
#endif

#include <iostream>

using namespace::std;

class LU_SOLVER
{
public:
  LU_SOLVER(): _pivots(NULL) { };
  ~LU_SOLVER()
  {
    if(_pivots){
      delete[] _pivots;
      _pivots = NULL;
    }
  }
  bool compute(const MATRIX& mat)
#if USING_OSX
  {
    _mat = mat;
    _rows = mat.rows();
    _cols = mat.cols();

    int smaller = (_rows < _cols) ? _rows : _cols;
    int info;
    
    if(_pivots){
      delete[] _pivots;
    }

    _pivots = new int[smaller];

    Real* address = &_mat(0, 0);
  
  #ifdef SINGLE_PRECISION
    sgetrf_((__CLPK_integer*)&_rows, (__CLPK_integer*)&_cols, address, (__CLPK_integer*)&_rows, (__CLPK_integer*)_pivots, (__CLPK_integer*)&info);
  #else
    dgetrf_((__CLPK_integer*)&_rows, (__CLPK_integer*)&_cols, address, (__CLPK_integer*)&_rows, (__CLPK_integer*)_pivots, (__CLPK_integer*)&info);
  #endif
    
    if (info < 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " LU factorization failed!" << endl;
      return false;
    }
    if (info > 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " _matrix is singular!" << endl;
      return false;
    }
    return true;
  }
#else
  {
    _mat.compute(mat);
  }
#endif

  void solve(const VECTOR& b, VECTOR& x)
#if USING_OSX
  {
    assert(_rows == b.size());
    x = b;

    if (_pivots == NULL)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " You need to call factorLU before you can call solveLU!" << endl;
      return;
    }

    // do we want to solve the transposed problem? (No)
    char transposed = 'N';

    // how many B's are we solving for? (one)
    int howManyBs = 1;

    int info;
    
    Real* address = &_mat(0, 0);

  #ifdef SINGLE_PRECISION
    sgetrs_(&transposed, (__CLPK_integer*)&_rows, (__CLPK_integer*)&howManyBs, address, (__CLPK_integer*)&_rows, (__CLPK_integer*)_pivots, &x[0], (__CLPK_integer*)&_rows, (__CLPK_integer*)&info);
  #else
    dgetrs_(&transposed, (__CLPK_integer*)&_rows, (__CLPK_integer*)&howManyBs, address, (__CLPK_integer*)&_rows, (__CLPK_integer*)_pivots, &x[0], (__CLPK_integer*)&_rows, (__CLPK_integer*)&info);
  #endif

    if (info != 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " LU solve failed!" << endl;
    }
  }
#else
  {
    x = _mat.solve(b);
  }
#endif

private:
  #if USING_OSX
  MATRIX _mat;
  #else
  Eigen::FullPivLU<MATRIX> _mat;
  #endif

  int _rows;
  int _cols;
  int* _pivots;
};

#endif
