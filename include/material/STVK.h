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
// STVK.h: interface for the STVK class.
//
//////////////////////////////////////////////////////////////////////

#ifndef STVK_H
#define STVK_H

#include <iostream>
#include <geometry/TET.h>
#include <material/MATERIAL.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// St-VK material class
//////////////////////////////////////////////////////////////////////
class STVK : public MATERIAL {

public:
  STVK(Real lambda, Real mu);
  ~STVK() {};

  MATRIX stiffness(TET& tet);
  MATRIX3 secondPiolaKirchhoff(const MATRIX3& F);
	Real strainEnergy(const MATRIX3& F);
  void stiffnessDensity(const Real* F, Real* stiffness);
  MATERIAL* copy();

private:
  Real _lambda;
  Real _mu;

  // work array
  MATRIX _pf_pF;

  void computeStresses(TET& tet);
};

#endif
