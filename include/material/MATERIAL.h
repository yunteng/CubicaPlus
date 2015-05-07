#ifndef MATERIAL_H
#define MATERIAL_H

#include <SETTINGS.h>
#include <geometry/TET.h>

using namespace std;

class TET;
//////////////////////////////////////////////////////////////////////
// Abstract material class
//////////////////////////////////////////////////////////////////////
class MATERIAL {

public:
  MATERIAL()
    : _materialName("")//, stiffnessDensityTime( 0.0 )
  {}
  virtual ~MATERIAL() {};

  // Set the parameter diagonal = true in cases when the deformation
  // gradient is explicitly known to be diagonal (cases in which
  // we diagonalize it ourselves).
  virtual MATRIX stiffness(TET& tet/*, bool diagonal = false*/) = 0;

  virtual Real strainEnergy(const MATRIX3& F) = 0;

  virtual MATRIX3 firstPiolaKirchhoff(const MATRIX3& F/*, bool diagonal = false*/)
  {
    return F * secondPiolaKirchhoff(F);
  }
  virtual MATRIX3 secondPiolaKirchhoff(const MATRIX3& F/*, bool diagonal = false*/)
  {
    return F.inverse() * firstPiolaKirchhoff(F);
  }
  virtual MATRIX3 cauchyStress(const MATRIX3& F/*, bool diagonal = false*/)
  {
    Real J = F.determinant();
    MATRIX3 S = secondPiolaKirchhoff(F);
    return (1.0 / J) * F * S * F.transpose();
  }

  virtual void stiffnessDensity(const Real* F, Real* stiffness/*, bool diagonal = false*/)
  {
    
  }

  // this is just a flattened version of first Piola-Kirchoff, and should
  // be renamed at some point
  //virtual void forceDensity(TET& tet, VEC3F* forces) = 0;
  // virtual void forceDensity(const Real* F, Real* forces/*, bool diagonal = false*/) = 0;

  // compute the invariant derivatives, 
  // \frac{\partial \Psi}{\partial I}
  // \frac{\partial \Psi}{\partial II}
  // \frac{\partial \Psi}{\partial III}
  // These are stacked in that order in the returned VEC3F
  VEC3F invariantPartials(const Real* F)
  {
    // The partials are solved for as follows:
    //
    // Using formula (6.23) from Bonet and Wood,
    // S = 2 \Psi_I I + 4 \Psi_II C + 2 J^2 \Psi_III C^-1
    // extract the first columns from S, 2I, 4C, and 2 J^2 C^-1
    // and stack \Psi_I, \Psi_II, \Psi_III into a vector.
    //
    // This is now a linear system that can be solved for all the partials.

    // put F into a matrix
    MATRIX3 F3x3;
    F3x3(0,0) = F[0];
    F3x3(1,0) = F[1];
    F3x3(2,0) = F[2];
    F3x3(0,1) = F[3];
    F3x3(1,1) = F[4];
    F3x3(2,1) = F[5];
    F3x3(0,2) = F[6];
    F3x3(1,2) = F[7];
    F3x3(2,2) = F[8];

    MATRIX3 S = secondPiolaKirchhoff(F3x3);
    MATRIX3 C = F3x3.transpose() * F3x3;
    MATRIX3 CInv = C.inverse();
    Real J = C.determinant();

    MATRIX3 A;
    A(0,0) = 2.0;
    A(1,0) = 0.0;
    A(2,0) = 0.0;
    A(0,1) = 4.0 * C(0,0);
    A(1,1) = 4.0 * C(1,0);
    A(2,1) = 4.0 * C(2,0);
    A(0,2) = 2.0 * J * J * CInv(0,0);
    A(1,2) = 2.0 * J * J * CInv(1,0);
    A(2,2) = 2.0 * J * J * CInv(2,0);
    VEC3F b;
    b[0] = S(0,0);
    b[1] = S(1,0);
    b[2] = S(2,0);
    VEC3F invariants = A.inverse() * b;
    return invariants;
  }

  void computePFPu(TET& tet, MATRIX& pFpu);

  // return a copy of this material object
  virtual MATERIAL* copy() = 0;
  string& materialName() { return _materialName; };


protected:
  string _materialName;
};

/*
Much of the code in this class and its subclasses was generated using Matlab and Maple.
For completeness, the generating code will be listed when possible. Some search-and-replace
is necessary to get the final code used in these classes, but the renamings necessary
should be self-evident.

The following are Matlab scripts used by all the material subclasses:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup_symbols_wrt_u.m
% setup symbols for when you want derivatives wrt u.
% This is more expensive than doing wrt just F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = sym('a'); b = sym('b'); c = sym('c'); d = sym('d');
e = sym('e'); f = sym('f'); g = sym('g'); h = sym('h');
i = sym('i'); j = sym('j'); k = sym('k'); l = sym('l');    
m = sym('m'); n = sym('n'); o = sym('o'); p = sym('p');
q = sym('q'); r = sym('r'); s = sym('s'); t = sym('t');
v = sym('u');    
n0 = [a b c];
n1 = [d e f];
n2 = [g h i];
n3 = [j k l];
u = mytrans([a b c d e f g h i j k l]);
x = mytrans([
    n1-n0
    n2-n0
    n3-n0
    ]);
Xinv = [
    m n o
    p q r
    s t v
    ];
F = x * Xinv;
C = mytrans(F)*F;
I_C = trace(C);
II_C = dblcon(C,C);
E = 0.5 * (C - eye(3));
J = det(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup_symbols_wrt_F.m
% setup symbols for when you only wants things wrt just F, not u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms F0 F1 F2 F3 F4 F5 F6 F7 F8;
F = [F0 F3 F6 
     F1 F4 F7 
     F2 F5 F8];
C = mytrans(F)*F;
I_C = trace(C);
II_C = dblcon(C,C);
E = 0.5 * (C - eye(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dblcon.m 
% symbolic double contraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = dblcon(A, B)
    [m,n] = size(A);
    s = 0;
    for i = 1:m
        for j = 1:n
            s = s + A(i,j)*B(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mytrans.m
% symbolic transpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = mytrans(M)
    [m,n] = size(M);
    N = sym('a')*eye(n,m);
    for i = 1:m
        for j = 1:n
            N(j,i) = M(i,j);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% codegen_density.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [forceDensity, stiffnessDensity] = codegen_density(strain_energy, F)
    Psi = strain_energy;
    F = reshape(F, 9, 1);
    display '-- computing stress force density --';
    Rd_sym = sym('DUMMY_TEMP')*ones(1,9);
    for i = 1:9
        Rd_sym(1,i) = - diff( Psi, F(i));
    end
    display '-- computing stiffness density --';
    Kd_sym = sym('DUMMY_TEMP')*ones(9*9,1);
    for i = 1:9
        for j = 1:9
            Kd_sym((i-1)*9+j) = diff( Rd_sym(1,i), F(j) );
        end
    end
    display '-- codegenning --';
    forceDensity = maple('codegen[C]', Rd_sym, 'optimized');
    stiffnessDensity = maple('codegen[C]', Kd_sym, 'optimized');
end
*/

#endif
