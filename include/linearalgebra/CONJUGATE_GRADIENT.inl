#include <util/SIMPLE_PARSER.h>
#include <util/MATRIX_UTIL.h>

template<class PRECONDITIONER, class INTEGRATOR>
CONJUGATE_GRADIENT<PRECONDITIONER,INTEGRATOR>::CONJUGATE_GRADIENT(INTEGRATOR* integrator, PRECONDITIONER* preconditioner):
  _integrator(integrator),
  _preconditioner(preconditioner),
  _totalFrames(0),
  _totalIterations(0)
{
  _maxIteration = SIMPLE_PARSER::getInt("cg max iterations", 10);
  _eps = SIMPLE_PARSER::getFloat("cg eps", 0.0001);
}
template<class PRECONDITIONER, class INTEGRATOR>
CONJUGATE_GRADIENT<PRECONDITIONER,INTEGRATOR>::~CONJUGATE_GRADIENT()
{

}
template<class PRECONDITIONER, class INTEGRATOR>
bool CONJUGATE_GRADIENT<PRECONDITIONER,INTEGRATOR>::solve(const VECTOR& b, VECTOR& x)
{
  if(b.size() == 0)
    return true;
  
  if(_preconditioner == NULL)
    return solveCG(b, x);
  else
    return solvePCG(b, x);
}
template<class PRECONDITIONER, class INTEGRATOR>
bool CONJUGATE_GRADIENT<PRECONDITIONER,INTEGRATOR>::solveCG(const VECTOR& b, VECTOR& x)
{
  x.conservativeResize(b.size());
  x.setZero();

  _integrator->getMatVecMult(x, _q);

  _direction = _residual = b - _q;

  Real deltaNew = _residual.squaredNorm();

  Real delta0 = deltaNew;

  bool converged = false;
  int i = 0;
  for(; i < _maxIteration && !converged; i++){

    _integrator->getMatVecMult(_direction, _q);

    Real alpha = _direction.dot(_q);
    if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;

    x += alpha * _direction;

    if(i > 0 && i % 25 == 0){
      _integrator->getMatVecMult(x, _residual);

      _residual *= -1;
      _residual += b;
    }
    else
      _residual -= alpha * _q;

    Real deltaOld = deltaNew;

    deltaNew = _residual.squaredNorm();

    Real beta = deltaNew / deltaOld;

    _direction *= beta;
    _direction += _residual;

    if(deltaNew <= _eps * delta0)
      converged = true;
  }

  if(converged){
    // cout << "Conjugate Gradient solver converged in " << i << " iterations" << endl;
    _totalIterations++;
    _totalFrames += i;
    _numberOfIterations.push_back(_totalFrames);
  }else
    cout << "Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;
  return converged;
}

template<class PRECONDITIONER, class INTEGRATOR>
bool CONJUGATE_GRADIENT<PRECONDITIONER,INTEGRATOR>::solvePCG(const VECTOR& b, VECTOR& x)
{
  x.conservativeResize(b.size());
  x.setZero();

  _integrator->getMatVecMult(x, _q);

  _residual = b - _q;

  _preconditioner->solve(_residual, _direction);

  Real deltaNew = _residual.dot(_direction);

  Real delta0 = deltaNew;

  bool converged = false;
  int i = 0;
  for(; i < _maxIteration && !converged; i++){
    _integrator->getMatVecMult(_direction, _q);

    Real alpha = _direction.dot(_q);
    if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;

    x += alpha * _direction;

    if(i > 0 && i % 25 == 0){
      _integrator->getMatVecMult(x, _residual);
      _residual *= -1;
      _residual += b;
    }
    else
      _residual -= alpha * _q;

    _preconditioner->solve(_residual, _s);

    Real deltaOld = deltaNew;

    deltaNew = _residual.dot(_s);


    Real beta = deltaNew / deltaOld;

    _direction *= beta;
    _direction += _s;

    if(deltaNew <= _eps * delta0)
      converged = true;
  }

  if(converged){
    _totalIterations++;
    _totalFrames += i;
    _numberOfIterations.push_back(i);
  }
  else
    cout << "Preconditioned Conjugate Gradient solver didn't converge in " << _maxIteration << " iterations, residual " << deltaNew << endl;

  return converged;
}
