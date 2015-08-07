/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <material/COROTATION.h>
#include <Eigen/SVD>
#include <util/MATRIX_UTIL.h>

//////////////////////////////////////////////////////////////////////
// Constructor for COROTATION
//////////////////////////////////////////////////////////////////////
COROTATION::COROTATION(Real lambda, Real mu) :
  _lambda(lambda), _mu(mu)
{
  _materialName.assign("COROTATION");
}

//////////////////////////////////////////////////////////////////////
// make a copy
//////////////////////////////////////////////////////////////////////
COROTATION* COROTATION::copy()
{
  COROTATION* material = new COROTATION(_lambda, _mu);
  return material;
}
// assume F is diagonal matrix
Real COROTATION::strainEnergy(const MATRIX3& F)
{
  MATRIX3 U, V, Fhat, R, S, L;
  decomposeF(F, U, Fhat, V, R, S, L);
  return strainEnergy(R, S);
}
// E = mu * ||F - R||_F + lambda / 2 * (S - I)
Real COROTATION::strainEnergy(const MATRIX3& R, const MATRIX3& S)
{
  MATRIX3 E = S - MATRIX3::Identity(); 
  
  Real energy = _mu * trace(E ^ E) + 0.5 * _lambda * (trace(E) * trace(E));
  
	return energy;
}
//////////////////////////////////////////////////////////////////////
// implementation of first PK stress tensor w.r.t. the singular values
//////////////////////////////////////////////////////////////////////
MATRIX3 COROTATION::firstPiolaKirchhoff(const MATRIX3& F/*, bool diagonal*/)
{
  MATRIX3 U, V, Fhat, R, S, L;
  decomposeF(F, U, Fhat, V, R, S, L);
  return firstPiolaKirchhoff(R, S);
}
MATRIX3 COROTATION::firstPiolaKirchhoff(const MATRIX3& R, const MATRIX3& S)
{
  MATRIX3 E = S - MATRIX3::Identity();
  MATRIX3 tmp = 2.0 * _mu * E;
  tmp = (tmp.array() + _lambda * trace(E)).matrix();
  // MATRIX3 firstPK = R * tmp;
  return R * tmp;
}

//////////////////////////////////////////////////////////////////////
// stiffness matrix implementation
//
// Does this need to be filtered? The tet being near inversion in the
// rest pose implies it is ill-conditioned, which shouldn't be happening
// with isosurface stuffing.
//////////////////////////////////////////////////////////////////////
MATRIX COROTATION::stiffness(TET& tet/*, bool diagonal*/)
{
  MATRIX3 F = tet.F();
    
  MATRIX3 U, V, Fhat, R, S, L;
  decomposeF(F, U, Fhat, V, R, S, L);
  
  // compute Bm matrix
  const VEC3F* b = tet.b();
  
  // get PFPu from the tet - only depends on rest state,
  // so don't need to update tet state at all
  MATRIX& pFpu = tet.PFPu();
  MATRIX stiffnessFinal(12,12);

  // compute each column of the matrix
  for (int y = 0; y < 12; y++)
  {
    // extract a column from pFpu (ie pretend 'x' is just a 
    // Kronecker delta)
    VECTOR deltaF(9);
    for (int z = 0; z < 9; z++)
      deltaF(z) = pFpu(z, y);
    
    MATRIX3 dF;
    MATRIX_UTIL::repackVec9(deltaF, dF);
    // MATRIX3 unrotatedDeltaP = material->firstPKDifferential(R, S, dF);
    MATRIX3 dFhat = R.transpose() * dF;
    MATRIX3 deltaP = rotatedFirstPKDifferential(L, dFhat); 
    
    deltaP = R * deltaP;
    
    VEC3F forceVecs[4];
    forceVecs[0] = deltaP * b[0];
    forceVecs[1] = deltaP * b[1];
    forceVecs[2] = deltaP * b[2];
    forceVecs[3] = deltaP * b[3];
   
    // copy result into stiffness column
    for (int z = 0; z < 4; z++)
      for (int a = 0; a < 3; a++)
        stiffnessFinal(z * 3 + a, y) = forceVecs[z][a];
  }

  return stiffnessFinal;
}
//////////////////////////////////////////////////////////////////////////////
// F = U * Fhat * VT
// R = U * VT
// S = V * Fhat * VT
//////////////////////////////////////////////////////////////////////////////
void COROTATION::decomposeF(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V, MATRIX3& R, MATRIX3& S, MATRIX3& L)
{
  svd(F, U, Fhat, V);

  R = U * V.transpose();
  S = V * Fhat * V.transpose();
  
  MATRIX3 Ld = (_mu * MATRIX3::Identity()) + (_lambda * trace(Fhat - MATRIX3::Identity()) - 2 * _mu) * ((trace(Fhat) * MATRIX3::Identity()) - Fhat).inverse();
  
  // cout << "Ld " << Ld << endl;
  // MATRIX3 t1 =  Ld - (_mu * MATRIX3::Identity());
  // cout << "Ld - mu * I " << t1 << endl;
  // Real coeff = (_lambda * trace(Fhat - MATRIX3::Identity()) - 2 * _mu);
  // t1 /= coeff;
  // MATRIX3 invt1 = t1.inverse();
  // cout << "(trace(Fhat) * MATRIX3::Identity()) - Fhat " << (trace(Fhat) * MATRIX3::Identity()) - Fhat << endl;
  // cout << "invt1 " << invt1 << endl;
  
  Ld(0, 0) = Ld(0, 0) > 0 ? Ld(0, 0) : 0;
  Ld(1, 1) = Ld(1, 1) > 0 ? Ld(1, 1) : 0;
  Ld(2, 2) = Ld(2, 2) > 0 ? Ld(2, 2) : 0;
  // for(int i = 0; i < 3; i++){
  //   if(Ld(i, i) < 0){
  //     Ld(i, i) = 0;
  //   }   
  // }

  L = V * Ld * V.transpose();
  
  
  // MATRIX3 newL = (_mu * MATRIX3::Identity()) + (_lambda * trace(S - MATRIX3::Identity()) - 2 * _mu) * ((trace(S) * MATRIX3::Identity()) - S).inverse();
  // cout << "newL - L " << newL - L << endl;
}
//////////////////////////////////////////////////////////////////////////////
// return the differential of the rotation matrix w.r.t dF
//////////////////////////////////////////////////////////////////////////////
MATRIX3 COROTATION::rotationDifferential(const MATRIX3& R, const MATRIX3& S, const MATRIX3& dF)
{
  // MATRIX3 W = R.transpose() * dF;
  MATRIX3 W = R ^ dF;
  VEC3F w(W(1, 2) - W(2, 1), W(2, 0) - W(0, 2), W(0, 1) - W(1, 0));
  // E = tr(S) * I - S
  MATRIX3 Einv = ((trace(S) * MATRIX3::Identity()) - S).inverse().eval();
  VEC3F r = Einv * w;
  return R * MATRIX_UTIL::cross(-r);
}
//////////////////////////////////////////////////////////////////////////////
// return the differential of the 1st PK stress tensor w.r.t. dF
//////////////////////////////////////////////////////////////////////////////
MATRIX3 COROTATION::firstPKDifferential(const MATRIX3& R, const MATRIX3& S, const MATRIX3& dF)
{
  MATRIX3 dR = rotationDifferential(R, S, dF);
  return (2.0 * _mu * dF) + (_lambda * trace(R.transpose() * dF) * R) + (_lambda * trace(S - MATRIX3::Identity()) - 2.0 * _mu) * dR;
}
MATRIX3 COROTATION::rotatedFirstPKDifferential(const MATRIX3& L, const MATRIX3& dFhat)
{
  MATRIX3 dFsym = 0.5 * (dFhat + dFhat.transpose());
  MATRIX3 dFskew = 0.5 * (dFhat - dFhat.transpose());
  
  MATRIX3 dPsym = (2 * _mu * dFsym) + (_lambda * trace(dFsym) * MATRIX3::Identity());

  VEC3F f(-dFskew(1, 2) + dFskew(2, 1), -dFskew(2, 0) + dFskew(0, 2), -dFskew(0, 1) + dFskew(1, 0));

  // MATRIX3 dPsym = _mu * (dFhat + dFhat.transpose());
  Real tr = _lambda * trace(dFhat);
  dPsym(0, 0) += tr;
  dPsym(1, 1) += tr;
  dPsym(2, 2) += tr;
  
  // VEC3F f(dFhat(2, 1) - dFhat(1, 2), dFhat(0, 2) - dFhat(2, 0), dFhat(1, 0) - dFhat(0, 1));
  
  // return dPsym + MATRIX3::cross(L * f);
  MATRIX3 dPskew = MATRIX_UTIL::cross(L * f);
  return dPsym + dPskew;
}

// void COROTATION::forceDensity(const Real* F, Real* forces/*, bool diagonal*/)
// { 
//   // MATRIX3 packedF = TET::repackF(F);
//   MATRIX3 dF;
//   MATRIX_UTIL::repackVec9(deltaF, dF);

//   MATRIX3 P = firstPiolaKirchhoff(packedF);
//   P *= -1.0;

//   forces[0] = P(0,0);
//   forces[1] = P(1,0);
//   forces[2] = P(2,0);
//   forces[3] = P(0,1);
//   forces[4] = P(1,1);
//   forces[5] = P(2,1);
//   forces[6] = P(0,2);
//   forces[7] = P(1,2);
//   forces[8] = P(2,2);
// }

//////////////////////////////////////////////////////////////////////////////
// Correct both Fhat and U if U contains a reflection
//////////////////////////////////////////////////////////////////////////////
void COROTATION::removeUReflection(MATRIX3& U, MATRIX3& Fhat)
{
  // find the smallest singular value
  int smallest = (Fhat(0,0) < Fhat(1,1)) ? 0 : 1;
  smallest = (Fhat(smallest, smallest) < Fhat(2,2)) ? smallest : 2;
 
  // negate it, and the corresponding column in U
  Fhat(smallest, smallest) *= -1.0;
  U(0, smallest) *= -1.0;
  U(1, smallest) *= -1.0;
  U(2, smallest) *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////
// Orthogonalize U if one of the singular values is near zero
//////////////////////////////////////////////////////////////////////////////
void COROTATION::orthogonalizeU(MATRIX3& U, MATRIX3& Fhat)
{
  // record the sign of U (ie whether it is a reflection) so that
  // the sign can be preserved.
  //
  // If it contains a reflection, it will be corrected in
  // removeUReflection, so make sure we don't double-handle it here.
  Real signU = (U.determinant() >= 0.0) ? 1.0 : -1.0;
  
  Real smallEigenvalue = 1e-4;

  // see if an eigenvalue is really small
  bool invalid[] = {false, false, false};
  int totalInvalid = 0;
  for (int x = 0; x < 3; x++)
    if (Fhat(x,x) < smallEigenvalue)
    {
      invalid[x] = true;
      totalInvalid++;
    }

  // if none are bad, do nothing
  if (totalInvalid == 0) return;

  // correct U according to the number of small eigenvalues
  if (totalInvalid == 1)
  {
    // figure out which column is bad
    int whichInvalid = 0;
    if (invalid[1]) whichInvalid = 1;
    if (invalid[2]) whichInvalid = 2;
    
    // get the other two columns
    VEC3F col0 = U.col((whichInvalid + 1) % 3);
    VEC3F col1 = U.col((whichInvalid + 2) % 3);

    // compute an orthogonal vector
    VEC3F col2 = col0.cross(col1);
    col2.normalize();

    // store it in U
    U(0, whichInvalid) = col2[0];
    U(1, whichInvalid) = col2[1];
    U(2, whichInvalid) = col2[2];

    // see if we changed any reflections
    if (signU * U.determinant() < 0.0)
    {
      U(0, whichInvalid) = -col2[0];
      U(1, whichInvalid) = -col2[1];
      U(2, whichInvalid) = -col2[2];
    }
    
    return;
  }
  else if (totalInvalid == 2)
  {
    // find the good vector
    int whichValid = 0;
    if (!invalid[1]) whichValid = 1;
    if (!invalid[2]) whichValid = 2;
    VEC3F validVec = U.col(whichValid);
   
    // find the smallest entry in the good vector
    int smallest = (fabs(validVec[0]) < fabs(validVec[1])) ? 0 : 1;
    smallest = (fabs(validVec[2]) < fabs(validVec[smallest])) ? 2 : smallest;
    
    // compute something arbitrarily orthogonal to
    // the good vector
    VEC3F secondValid;
    int next = (smallest + 1) % 3;
    int nextnext = (smallest + 2) % 3;
    secondValid[smallest] = 0.0;
    secondValid[next] = -validVec[nextnext];
    secondValid[nextnext] = validVec[next];

    // compute a third good vector
    VEC3F thirdValid = validVec ^ secondValid;
    thirdValid.normalize();

    // copy the good vectors into U
    next = (whichValid + 1) % 3;
    nextnext = (whichValid + 2) % 3;
    U(0, next) = secondValid[0];
    U(1, next) = secondValid[1];
    U(2, next) = secondValid[2];
    U(0, nextnext) = thirdValid[0];
    U(1, nextnext) = thirdValid[1];
    U(2, nextnext) = thirdValid[2];

    // see if we changed any reflections
    if (signU * U.determinant() < 0.0)
    {
      // negate the one with the smallest F, which is the one
      // removeUReflection will be negating again
      if (Fhat(next,next) < Fhat(nextnext, nextnext))
      {
        U(0, next) = -secondValid[0];
        U(1, next) = -secondValid[1];
        U(2, next) = -secondValid[2];
      }
      else
      {
        U(0, nextnext) = -thirdValid[0];
        U(1, nextnext) = -thirdValid[1];
        U(2, nextnext) = -thirdValid[2];
      }
    }
    return;
  }

  // all three are bad, just use identity
  U = MATRIX3::Identity();
}

//////////////////////////////////////////////////////////////////////////////
// Get the SVD of F
//////////////////////////////////////////////////////////////////////////////
void COROTATION::svd(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V)
{
  // compute the SVD
  MATRIX3 Fnormal3 = F.transpose() * F;
  MATRIX Fnormal(Fnormal3);

  // Eigen::JacobiSVD<MATRIX3> eigenSystem(Fnormal3, Eigen::ComputeFullV);
  // VEC3F eigenvalues = eigenSystem.singularValues();

  VEC3F eigenvalues;

  MATRIX_UTIL::eigensystem3x3(Fnormal3, eigenvalues, V);
  
  // populate the diagonal
  Fhat.setZero();
  Fhat(0,0) = sqrt(eigenvalues(0));
  Fhat(1,1) = sqrt(eigenvalues(1));
  Fhat(2,2) = sqrt(eigenvalues(2));

#ifdef _WIN32  
  if (_isnan(Fhat(0,0))) Fhat(0,0) = 0.0;
  if (_isnan(Fhat(1,1))) Fhat(1,1) = 0.0;
  if (_isnan(Fhat(2,2))) Fhat(2,2) = 0.0;
#else
  if (std::isnan(Fhat(0,0))) Fhat(0,0) = 0.0;
  if (std::isnan(Fhat(1,1))) Fhat(1,1) = 0.0;
  if (std::isnan(Fhat(2,2))) Fhat(2,2) = 0.0;
#endif

  // populate V
  // eigenvectors.copiesInto(V);
  // V = eigenSystem.matrixV();
  
  // if V is a reflection, multiply a column by -1
  // to ensure it is a rotation
  Real detV = V.determinant();
  if (detV < 0.0)
  {
    V(0,0) *= -1.0;
    V(1,0) *= -1.0;
    V(2,0) *= -1.0;
  }

  // compute U
  U.setZero();
  U(0,0) = (Fhat(0,0) > 0.0f) ? 1.0f / Fhat(0,0) : 0.0f;
  U(1,1) = (Fhat(1,1) > 0.0f) ? 1.0f / Fhat(1,1) : 0.0f;
  U(2,2) = (Fhat(2,2) > 0.0f) ? 1.0f / Fhat(2,2) : 0.0f;
  U = F * V * U;
  orthogonalizeU(U, Fhat);
 
  // correct any reflections in U to ensure it is a rotation
  if (F.determinant() < 0.0)
    removeUReflection(U, Fhat);

  
  ///////////////////////////////////////////////////////////////
  // All error checking starting here. Fast version would macro
  // this out
  ///////////////////////////////////////////////////////////////
  
  // basic error checking
  /*
  if (fabs(U.determinant() - 1.0) > 1e-4 || 
      fabs(V.determinant() - 1.0) > 1e-4)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " U and V of SVD have drifted from pure rotation!" << endl;
    cout << " det(U): " << U.determinant() << endl;
    cout << " det(V): " << V.determinant() << endl;
  }
  
  // sanity check that this decomposition does hand back F when combined
  MATRIX3 sanity = U * Fhat * V.transpose();
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++)
    {
      Real diff = sanity(x,y) - F(x,y);
      // if (F(x,y) > 0.0)
        // diff = diff / F(x,y);
      if (diff > 1e-4)
      {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
        cout << "SVD and original matrix differ!" << endl;
        cout << " relative diff: " << diff << " sanity: " << sanity(x,y) << " original: " << F(x,y) << endl;
        cout << " SVD product: " << sanity << endl;
        cout << " original: " << F << endl;
        cout << " normal F: " << Fnormal3 << endl;
        cout << " eigenvalues: " << eigenvalues << endl;
        cout << " eigenvectors: " << eigenSystem.matrixV() << endl;
        cout << " eigenproduct: " << eigenSystem.matrixV() * eigenvalues.asDiagonal() * eigenSystem.matrixV().transpose();
        // exit(0);
      }
    }
  */
}
