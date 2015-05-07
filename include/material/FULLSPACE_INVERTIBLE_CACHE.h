#ifndef FULLSPACE_INVERTIBLE_CACHE_H
#define FULLSPACE_INVERTIBLE_CACHE_H

#include <geometry/TET_MESH.h>
#include <material/INVERTIBLE.h>

class FULLSPACE_INVERTIBLE_CACHE
{
public:
  FULLSPACE_INVERTIBLE_CACHE(TET_MESH* tetMesh);
  ~FULLSPACE_INVERTIBLE_CACHE();

  void cacheDecompositions();
  Real computeElasticEnergy();
  VECTOR& computeInternalForce();
  COO_MATRIX& computeStiffnessMatrix();
  MATRIX reduceStiffness(MATRIX& U)
  {
    return U.transpose() * _tetMesh->stiffnessMatrix().rightMult(U);
  }

  void computeIndividualTetInternalForces(VECTOR& output);
  void computeIndividualTetInternalForces(vector<int>& tetIDs, VECTOR& output);

  // dummy function to keep this class compatible with
  // the subspace version
  MATRIX computeReducedStiffnessMatrix()
  {
    cout << __FILE__ << " " << __FUNCTION__ << endl
         << " You shouldn't be calling this!!!" << endl;
    exit(0);
    return MATRIX();
  }
  void cacheKeyTetTransforms(vector<MATRIX3>& dummy)
  {
    cout << __FILE__ << " " << __FUNCTION__ << endl
         << " You shouldn't be calling this!!!" << endl;
    exit(0);
  }

private:
  TET_MESH* _tetMesh;
  int _totalCores;
  vector<MATRIX3> _Us;
  vector<MATRIX3> _Vs;
  vector<MATRIX3> _Fhats;
  vector<MATRIX9> _stiffnesses;

  VECTOR* _RCopies;
};
#endif
