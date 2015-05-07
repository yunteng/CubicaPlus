#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include <SETTINGS.h>
#include <iostream>
#include <util/TIMING_BREAKDOWN.h>
#include <util/SIMPLE_PARSER.h>
#include <util/MATRIX_UTIL.h>
#include <linearalgebra/COO_MATRIX.h>

using namespace::std;

template<class INTEGRATOR, class SOLVER>
class NEWTON_RAPHSON
{
protected:
  INTEGRATOR* _integrator;
  SOLVER*     _solver;
  bool        _constHessian;
  bool        _constCollision;
  Real        _energy;
  Real        _minEnergy;
  VECTOR      _bestSolution;
  VECTOR      _solution;
  VECTOR      _delta;
  int         _maxNewtonIterations;
  Real        _eps;

public:
  NEWTON_RAPHSON(INTEGRATOR* integrator, SOLVER* solver):
    _integrator(integrator),
    _solver(solver)
  {
    _constHessian = SIMPLE_PARSER::getBool("constant hessian", false);
    _constCollision = SIMPLE_PARSER::getBool("constant collision jacobian", true);
    _maxNewtonIterations = SIMPLE_PARSER::getInt("max newton iterations", 75);
    _eps = SIMPLE_PARSER::getFloat("newton eps", 0.01);
    _delta.resize(_integrator->getPosition().size());
    _delta.setZero();
  }
  ~NEWTON_RAPHSON()
  {
  }
  bool run()
  {
    // TIMING_BREAKDOWN::tic();
    _integrator->initializeImplicitStep();
    // TIMING_BREAKDOWN::toc("Initialize Implicit Step");

    int iteration = 0;
    bool converged = false;

    Real currentEnergy = 0;
    Real currentResidual = 0;
    Real previousEnergy = 0;
    Real previousResidual = 0;

    bool redo = false;

    _solution = _integrator->getPosition();

    for(; iteration < _maxNewtonIterations; iteration++){
      computeLinearSystem(iteration == 0 || !_constHessian || redo, !_constCollision);

      currentResidual = _integrator->gradient().norm();
      currentEnergy = _integrator->energy();

      if(iteration == 0){
        previousEnergy = _integrator->energy();
        previousResidual = currentResidual;
      }
      // a hacky way to detect whether the Newton process has diverged
      else if(currentEnergy >= 4 * previousEnergy || currentResidual >= 4 * previousResidual){
        if(redo || !_constHessian){
          cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
          cout << "Newton iteration diverged! Terminate early..." << endl;
          return false;
        }else{
          cout << "Newton iteration diverged! redo with non-constant Hessain!!!" << endl;
          redo = true;
          iteration = -1;
          _integrator->resetState();
          continue;
        }
      }
      previousEnergy = currentEnergy;
      previousResidual = currentResidual;

      if(iteration == 0 || _minEnergy > currentEnergy){
        _minEnergy = currentEnergy;
        _bestSolution = _solution;
      }
      
      if(currentResidual < _eps){
        converged = true;
        break;
      }

      TIMING_BREAKDOWN::tic();
      _solver->solve(_integrator->gradient(), _delta);
      TIMING_BREAKDOWN::toc("Solve Linear System");

      _solution -= _delta;

      TIMING_BREAKDOWN::tic();
      _integrator->setPosition(_solution);
      TIMING_BREAKDOWN::toc("Set Position"); 
    }
    if(currentEnergy > _minEnergy){
      TIMING_BREAKDOWN::tic();
      _integrator->setPosition(_bestSolution);
      TIMING_BREAKDOWN::toc("Set Position");
    }
    
    TIMING_BREAKDOWN::tic();
    _integrator->finalizeImplicitStep();
    TIMING_BREAKDOWN::toc("Finalize Implicit Step");

    return converged;
  }
  
  void computeLinearSystem(bool recomputeHessian, bool recomputeCollision)
  {
    TIMING_BREAKDOWN::tic();
    _integrator->computeMaterialCache();
    TIMING_BREAKDOWN::toc("Compute Material Cache");

    TIMING_BREAKDOWN::tic();
    _integrator->computeSystemEnergy();
    TIMING_BREAKDOWN::toc("Compute System Energy");

    TIMING_BREAKDOWN::tic();
    if(recomputeHessian){
      _integrator->computeSystemMatrix();
    }
    else if(recomputeCollision){
      _integrator->computeCollisionMatrices();
    }

    TIMING_BREAKDOWN::toc("Compute System Matrix");
    
    TIMING_BREAKDOWN::tic();
    _integrator->computeSystemForce();
    TIMING_BREAKDOWN::toc("Compute System Force");

    _energy = _integrator->energy();

    cout << "energy: " << _integrator->energy() << " residual norm: " << _integrator->gradient().norm() << endl;
  }
};

#endif
