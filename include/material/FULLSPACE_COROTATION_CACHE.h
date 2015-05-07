#ifndef FULLSPACE_COROTATION_CACHE_H
#define FULLSPACE_COROTATION_CACHE_H

#include <geometry/TET_MESH.h>
#include <material/COROTATION.h>

class FULLSPACE_COROTATION_CACHE
{
public:
  FULLSPACE_COROTATION_CACHE(TET_MESH* tetMesh);
  ~FULLSPACE_COROTATION_CACHE();

  void cacheDecompositions(bool surfaceOnly = false);
  Real computeElasticEnergy();
  VECTOR& computeInternalForce();
  VECTOR& computeSurfaceInternalForce();
  COO_MATRIX& computeStiffnessMatrix();
  BLOCK_SPARSE_MATRIX& computeSurfaceStiffnessMatrix();

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

  vector<MATRIX3> _Rs;
  vector<MATRIX3> _Ss;
  vector<MATRIX3> _Ls;

  VECTOR* _RCopies;
};
#endif
