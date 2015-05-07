#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <SETTINGS.h>
#include <iostream>

using namespace::std;

template<class PRECONDITIONER, class INTEGRATOR>
class CONJUGATE_GRADIENT
{
public:
  CONJUGATE_GRADIENT(INTEGRATOR* integrator, PRECONDITIONER* preconditioner = NULL);
  ~CONJUGATE_GRADIENT();

  int& maxIteration() { return _maxIteration; };
  Real& eps()   { return eps; };

  bool solve(const VECTOR& b, VECTOR& x);
  bool solveCG(const VECTOR& b, VECTOR& x);
  bool solvePCG(const VECTOR& b, VECTOR& x);

  vector<int>& numberOfIterations() { return _numberOfIterations; };

protected:
  INTEGRATOR* _integrator;
  PRECONDITIONER* _preconditioner;

  int _maxIteration;
  Real _eps;

  int _totalFrames;
  int _totalIterations;
  vector<int> _numberOfIterations;

  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  VECTOR _s;
};

#include "CONJUGATE_GRADIENT.inl"
#endif
