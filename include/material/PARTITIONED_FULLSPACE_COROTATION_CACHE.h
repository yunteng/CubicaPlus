/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef PARTITIONED_FULLSPACE_COROTATION_CACHE_H
#define PARTITIONED_FULLSPACE_COROTATION_CACHE_H

#include <geometry/TET_MESH.h>
#include <material/COROTATION.h>

class PARTITIONED_FULLSPACE_COROTATION_CACHE
{
public:
  PARTITIONED_FULLSPACE_COROTATION_CACHE(TET_MESH* tetMesh);
  ~PARTITIONED_FULLSPACE_COROTATION_CACHE();

  void cacheDecompositions();
  void cachePartialDecompositions();

  Real computeElasticEnergy();
  
  VECTOR& computeInternalForce();

  BLOCK_COO_MATRIX& computeStiffnessMatrix();

  void computePartialStiffnessMatrix(int partition, COO_MATRIX& diagMat, COO_MATRIX& offDiagMat, COO_MATRIX& diagMat2);

  VECTOR& computePartialInternalForce(int partition);
  vector<VECTOR>& computePartialInternalForce();

  vector<VECTOR>& internalForces() { return _internalForces; };

private:
  TET_MESH* _tetMesh;

  vector<VECTOR> _internalForces;
  vector<MATRIX3> _Rs;
  vector<MATRIX3> _Ss;
  vector<MATRIX3> _Ls;
};
#endif
