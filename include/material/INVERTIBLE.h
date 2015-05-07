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
// INVERTIBLE.h: interface for the INVERTIBLE class.
//
//////////////////////////////////////////////////////////////////////

#ifndef INVERTIBLE_H
#define INVERTIBLE_H

#include <iostream>
#include <material/MATERIAL.h>
#include <util/TIMER.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Implements invertible finite elements from [Irving et al. 2004]
//
// The approach here is to filter the deformation gradient as it
// comes and pass it with clamped singular values to a concrete 
// MATERIAL subclass for evaluation.
//////////////////////////////////////////////////////////////////////
class INVERTIBLE : public MATERIAL {

public:
  static int WHICHCASE;
  static int NUMINVERTED;

  //INVERTIBLE(MATERIAL* material, Real epsilon = 0.1);
  INVERTIBLE(MATERIAL* material, bool doFind = true);
  ~INVERTIBLE() { delete _material; };

  // return a copy of this material object
  MATERIAL* copy();
  Real& epsilon() { return _epsilon; };

  // return the material being clamped
  MATERIAL* material() { return _material; };

	Real strainEnergy(const MATRIX3& F);

  MATRIX stiffness(TET& tet);
  MATRIX3 firstPiolaKirchhoff(const MATRIX3& F);
  MATRIX3 firstPiolaKirchhoff(const MATRIX3& U, const MATRIX3& Fhat, const MATRIX3& V);
  MATRIX3 secondPiolaKirchhoff(const MATRIX3& F);

  void stiffnessDensity(const Real* F, Real* stiffness);

  // diagonalized version
  void stiffnessDensity(const MATRIX3& U, 
                        const MATRIX3& Fhat,
                        const MATRIX3& V,
                        Real* stiffness);

  // diagonalize the deformation gradient
  void diagonalizeF(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);


  // access counter of how many tets were inverted
  static int& inversions() { return _inversions; };

  Real &stiffnessDensityTime() { return _stiffnessDensityTime; }
  Real &materialStiffnessDensityTime() { return _materialStiffnessDensityTime; }

private:
  MATERIAL* _material;
  Real _epsilon;  // the threshold at which to start clamping F
  Real _implicitEpsilon;  // the threshold at which to start clamping the jacobian
 
  // functions to take the SVD of deformation gradient F
  void removeUReflection(MATRIX3& U, MATRIX3& Fhat);
  void orthogonalizeU(MATRIX3& U, MATRIX3& Fhat);
  void svd(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);

  Real findRoot(Real find);

  void findClamping(bool scanRightLeft = false);

  Real _slope;
  Real _intercept;
  static int _inversions;

  Real _stiffnessDensityTime;
  Real _materialStiffnessDensityTime;
};

#endif


