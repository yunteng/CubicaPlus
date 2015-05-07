/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// COROTATION.h: interface for the COROTATION class.
//
//////////////////////////////////////////////////////////////////////

#ifndef COROTATION_H
#define COROTATION_H

#include <iostream>
#include <material/MATERIAL.h>
#include <geometry/TET.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// St-VK material class
//////////////////////////////////////////////////////////////////////
class COROTATION : public MATERIAL {

public:
  COROTATION(Real lambda, Real mu);
  ~COROTATION() {};

  Real lambda() { return _lambda; }
  Real mu() { return _mu; }

  MATRIX stiffness(TET& tet/*, bool diagonal = false*/);

  COROTATION* copy();

  Real strainEnergy(const MATRIX3& F);
  Real strainEnergy(const MATRIX3& R, const MATRIX3& S);
  
  MATRIX3 firstPiolaKirchhoff(const MATRIX3& F/*, bool diagonal = true*/);
  MATRIX3 firstPiolaKirchhoff(const MATRIX3& R, const MATRIX3& S);
  
  MATRIX3 firstPKDifferential(const MATRIX3& R, const MATRIX3& S, const MATRIX3& dF);
  MATRIX3 rotatedFirstPKDifferential(const MATRIX3& L, const MATRIX3& dFhat);
  
  // void forceDensity(TET& tet, VEC3F* forces/*, bool diagonal = false*/)
    // { cout << __FILE__ << " " << __LINE__ << " : UNIMPLEMENTED " << endl; };
  // void forceDensity(const Real* F, Real* forces/*, bool diagonal = false*/);

  // void stiffnessDensity(const Real* F, Real* stiffness/*, bool diagonal = true*/);
  void decomposeF(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V, MATRIX3& R, MATRIX3& S, MATRIX3& L);
  
private:
  Real _lambda;
  Real _mu;
  
  void svd(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);
  void removeUReflection(MATRIX3& U, MATRIX3& Fhat);
  void orthogonalizeU(MATRIX3& U, MATRIX3& Fhat);
  
  MATRIX3 rotationDifferential(const MATRIX3& R, const MATRIX3& S, const MATRIX3& dF);
};

#endif
