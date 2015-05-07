#ifndef SUBSPACE_COROTATION_CACHE_H
#define SUBSPACE_COROTATION_CACHE_H

#include <geometry/SUBSPACE_TET_MESH.h>
#include <material/COROTATION.h>

class SUBSPACE_COROTATION_CACHE
{
public:
  SUBSPACE_COROTATION_CACHE(SUBSPACE_TET_MESH* tetMesh);
  ~SUBSPACE_COROTATION_CACHE();

  void cacheDecompositions();
  void cacheKeyTetTransforms(vector<MATRIX3>& transformMatrix);
  Real computeElasticEnergy();
  VECTOR& computeInternalForce();
  MATRIX& computeReducedStiffnessMatrix();

  COO_MATRIX computeStiffnessMatrix()
  {
    cout << __FILE__ << " " << __FUNCTION__ << endl
         << " You shouldn't be calling this!!!" << endl;
    exit(0);
    return COO_MATRIX();
  }
  // dummy function to keep this class compatible with
  // the fullspace version
  MATRIX reduceStiffness(MATRIX& U)
  {
    cout << __FILE__ << " " << __FUNCTION__ << endl
         << " You shouldn't be calling this!!!" << endl;
    exit(0);
    return MATRIX();
  }

private:
  SUBSPACE_TET_MESH* _tetMesh;

  vector<MATRIX> _keyTetUs;
  vector<SpMat> _keyTetTransforms;

  vector<MATRIX3> _Rs;
  vector<MATRIX3> _Ss;
  vector<MATRIX3> _Ls;

  VECTOR* _RCopies;
  MATRIX* _KCopies;
};
#endif
