#ifndef HYBRID_PRECONDITIONER_H
#define HYBRID_PRECONDITIONER_H
#include <SETTINGS.h>
#include <iostream>
#include <linearalgebra/LU_SOLVER.h>
#if USING_OPENMP
#include <omp.h>
#endif

using namespace::std;
template<class INTEGRATOR>
class HYBRID_PRECONDITIONER
{
public:
  HYBRID_PRECONDITIONER(INTEGRATOR* integrator):
    _integrator(integrator),
    _diag(integrator->hessianDiagonal()), _dense(integrator->reducedKiiInv())
  {

  }
  ~HYBRID_PRECONDITIONER() {};

  void solve(const VECTOR& b, VECTOR& x)
  {
    x.conservativeResize(b.size());
    
    x.head(_diag.size()) = _diag.asDiagonal().inverse() * b.head(_diag.size());

    int leftOver = b.size() - _diag.size();

    if(_integrator->hasEntireFullPartitions()){
      _tailb = b.tail(leftOver);
      _integrator->removeZeroVectors(_tailb, _prunedTailb);
      _dense.solve(_prunedTailb, _prunedTailx);
      _integrator->fillInZeroVectors(_prunedTailx, _tailx);
      x.tail(leftOver) = _tailx;
    }else{
      _tailb = b.tail(leftOver);
      _dense.solve(_tailb, _tailx);
      x.tail(leftOver) = _tailx;
    }
  }

private:
  INTEGRATOR* _integrator;
  const VECTOR& _diag;
  LU_SOLVER& _dense;

  VECTOR _tailb;
  VECTOR _prunedTailb;
  VECTOR _prunedTailx;
  VECTOR _tailx;
};
#endif
