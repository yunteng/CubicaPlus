/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef JACOBI_PRECONDITIONER_H
#define JACOBI_PRECONDITIONER_H
#include <SETTINGS.h>
#include <iostream>
#if USING_OPENMP
#include <omp.h>
#endif

using namespace::std;
class JACOBI_PRECONDITIONER
{
public:
  JACOBI_PRECONDITIONER(VECTOR& diagonal):
    _diag(diagonal)
  {

  }
  ~JACOBI_PRECONDITIONER() {};

  void solve(const VECTOR& b, VECTOR& x)
  {
    x = (b.array() / _diag.array()).matrix();
  }

private:
  const VECTOR& _diag;
};
#endif
