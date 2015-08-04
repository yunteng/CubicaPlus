/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef PARTITIONED_SUBSPACE_COROTATION_CACHE_H
#define PARTITIONED_SUBSPACE_COROTATION_CACHE_H

#include <geometry/SUBSPACE_TET_MESH.h>
#include <material/COROTATION.h>

class PARTITIONED_SUBSPACE_COROTATION_CACHE
{
public:
  PARTITIONED_SUBSPACE_COROTATION_CACHE(SUBSPACE_TET_MESH* tetMesh);
  ~PARTITIONED_SUBSPACE_COROTATION_CACHE();

  void setCoupledWithFullspace(bool val) {_coupledWithFullspace = true; }

  void cacheDecompositions();
  void cacheKeyTetTransforms(vector<MATRIX3>& transformMatrix);
  Real computeElasticEnergy();
  VECTOR& computeInternalForce();
  VECTOR& computeUnprojectedInternalForce();
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

  void registerFullsimTets();

private:
  SUBSPACE_TET_MESH* _tetMesh;

  bool _coupledWithFullspace;

  vector<vector<MATRIX> > _keyTetUs;
  
  vector<SpMat> _keyTetTransforms;

  vector<MATRIX3> _Rs;
  vector<MATRIX3> _Ss;
  vector<MATRIX3> _Ls;

  VECTOR _isFullsimTets;
};
#endif
