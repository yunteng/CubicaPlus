/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <material/SUBSPACE_COROTATION_CACHE.h>
#include <util/MATRIX_UTIL.h>
#include <util/TIMING_BREAKDOWN.h>

#if USING_SUBSPACE_OPENMP
#include <omp.h>
#endif

SUBSPACE_COROTATION_CACHE::SUBSPACE_COROTATION_CACHE(SUBSPACE_TET_MESH* tetMesh):
  _tetMesh(tetMesh)
{
  int totalKeyTets = _tetMesh->keyTets().size();
  _Rs.resize(totalKeyTets);
  _Ss.resize(totalKeyTets);
  _Ls.resize(totalKeyTets);
  _keyTetTransforms.resize(totalKeyTets);

  for(int x = 0; x < totalKeyTets; x++){
    _keyTetTransforms[x].resize(12, 12);
    _keyTetTransforms[x].reserve(Eigen::VectorXi::Constant(12, 3));
  }
  _RCopies = new VECTOR[_tetMesh->totalCores()];
  _KCopies = new MATRIX[_tetMesh->totalCores()];
  for(int x = 0; x < _tetMesh->totalCores(); x++){
    _RCopies[x].resize(_tetMesh->rank());
    _KCopies[x].resize(_tetMesh->rank(), _tetMesh->rank());
  }
  _keyTetUs.clear();
  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();

  for(int x = 0; x < totalKeyTets; x++){
    MATRIX tetBasis = _tetMesh->tetBasis(keyTetsIDs[x]);
    TET& tet = tets[keyTetsIDs[x]];
    bool containsSurface = false;

    _keyTetUs.push_back(tetBasis);
  }
}
SUBSPACE_COROTATION_CACHE::~SUBSPACE_COROTATION_CACHE()
{
  delete[] _RCopies;
  delete[] _KCopies;
}
void SUBSPACE_COROTATION_CACHE::cacheDecompositions()
{
  TIMING_BREAKDOWN::tic();

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();
  _tetMesh->updateFullMesh();
  _tetMesh->generateKeyTetsF();
  vector<MATRIX3>& Fs = _tetMesh->F();

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_SUBSPACE_OPENMP
    const int id = omp_get_thread_num();
    #else
    const int id = 0;
    #endif

    #if USING_SUBSPACE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int x = 0; x < keyTetsIDs.size(); x++){
      int materialIndex = tets[keyTetsIDs[x]].materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
      MATRIX3& F = Fs[x];
      MATRIX3 U, V, Fhat;
      material->decomposeF(F, U, Fhat, V, _Rs[x], _Ss[x], _Ls[x]);
    }
  }
  TIMING_BREAKDOWN::toc("Cache Subspace Material Decompositions");
  
}

void SUBSPACE_COROTATION_CACHE::cacheKeyTetTransforms(vector<MATRIX3>& transformMatrix)
{
  TIMING_BREAKDOWN::tic();

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();

#if USING_SUBSPACE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int x = 0; x < keyTetsIDs.size(); x++){
    TET& tet = tets[keyTetsIDs[x]];

    for(int y = 0; y < 4; y++){
      int vertexID = _tetMesh->vertexID(tet.vertices[y]);
      if(!_tetMesh->isConstrained(vertexID)){
        for(int i = 0; i < 3; i++)
          for(int j = 0; j < 3; j++){
             _keyTetTransforms[x].coeffRef(y * 3 + i, y * 3 + j) = transformMatrix[vertexID](i, j);
          }
      }
    }
  }

  TIMING_BREAKDOWN::toc("Cache Key Tet Transforms");
  
}


VECTOR& SUBSPACE_COROTATION_CACHE::computeInternalForce()
{
  TIMING_BREAKDOWN::tic();

  VECTOR& R = _tetMesh->reducedInternalForce();
  if(R.size() != _tetMesh->rank())
    R.resize(_tetMesh->rank());
  
  R.setZero();

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();
  vector<Real>& keyTetWeights = _tetMesh->keyWeights();

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_SUBSPACE_OPENMP
    const int id  = omp_get_thread_num();

#else
    const int id  = 0;
#endif

    _RCopies[id].setZero();

  #if USING_SUBSPACE_OPENMP
  #pragma omp for schedule(static)
  #endif
  // populate the forces
    for (unsigned int x = 0; x < keyTetsIDs.size(); x++)
    {
      int tetID = keyTetsIDs[x];
      // compute the forces
      VECTOR forces(12);
      const VEC3F* b = tets[tetID].b();

      int materialIndex = tets[tetID].materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[x], _Ss[x]);
      forces.segment<3>(0) = firstPK * b[0];
      forces.segment<3>(3) = firstPK * b[1];
      forces.segment<3>(6) = firstPK * b[2];
      forces.segment<3>(9) = firstPK * b[3];

      if(_tetMesh->isSkinningBasis()){
        _RCopies[id] += keyTetWeights[x] * _keyTetUs[x].transpose() * _keyTetTransforms[x].transpose() * forces;
      }else{
        _RCopies[id] += keyTetWeights[x] * _keyTetUs[x].transpose() * forces;
      }
    }
  }
  // merge all the copies
  for (int x = 0; x < _tetMesh->totalCores(); x++)
    R += _RCopies[x];

  TIMING_BREAKDOWN::toc("Compute Reduced Internal Force");
  return R;
}

Real SUBSPACE_COROTATION_CACHE::computeElasticEnergy()
{
  TIMING_BREAKDOWN::tic();

  Real elasticEnergy = 0.0;

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();
  vector<Real>& keyTetWeights = _tetMesh->keyWeights();

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_OPENMP
    const int id = omp_get_thread_num();
    #else
    const int id = 0;
    #endif

    #if USING_SUBSPACE_OPENMP
    #pragma omp for schedule(static) reduction(+ : elasticEnergy)
    #endif
    for (int x = 0; x < keyTetsIDs.size(); x++)
    {
      int tetID = keyTetsIDs[x];

      int materialIndex = tets[tetID].materialIndex();
      
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      elasticEnergy += keyTetWeights[x] * tets[tetID].restVolume() * material->strainEnergy(_Rs[x], _Ss[x]);
    }
  }
  TIMING_BREAKDOWN::toc("Compute Reduced Elastic Energy");
  return elasticEnergy;
}

MATRIX& SUBSPACE_COROTATION_CACHE::computeReducedStiffnessMatrix()
{
  TIMING_BREAKDOWN::tic();

  MATRIX& stiffness = _tetMesh->reducedStiffnessMatrix();

  if(stiffness.rows() != _tetMesh->rank())
    stiffness.resize(_tetMesh->rank(), _tetMesh->rank());
  stiffness.setZero();

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();
  vector<Real>& keyTetWeights = _tetMesh->keyWeights();


#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  { 
    #if USING_SUBSPACE_OPENMP
    const int id  = omp_get_thread_num();
    #else
    const int id = 0;
    #endif
    _KCopies[id].setZero();
    
    #if USING_SUBSPACE_OPENMP
    #pragma omp for  schedule(static)
    #endif
    for (unsigned int x = 0; x < keyTetsIDs.size(); x++)
    {
      MATRIX tetStiffness(12, 12);
      tetStiffness.setZero();

      int tetID = keyTetsIDs[x];

      TET& tet = tets[tetID];
      int materialIndex = tet.materialIndex();
      // compute Bm matrix
      const VEC3F* b = tet.b();

      int indices[] = {_tetMesh->vertexID(tet.vertices[0]),
                       _tetMesh->vertexID(tet.vertices[1]),
                       _tetMesh->vertexID(tet.vertices[2]),
                       _tetMesh->vertexID(tet.vertices[3])};
      
      // get PFPu from the tet - only depends on rest state,
      // so don't need to update tet state at all
      MATRIX& pFpu = tet.PFPu();
      
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
      
      MATRIX3& R = _Rs[x];
      MATRIX3& S = _Ss[x];
      MATRIX3& L = _Ls[x];
      // compute each column of the matrix
      for (int y = 0; y < 12; y++)
      {
        int col = indices[y / 3] * 3 + y % 3;
        if(col >= _tetMesh->dofs())
          continue;

        // extract a column from pFpu (ie pretend 'x' is just a 
        // Kronecker delta)
        VECTOR deltaF(9);
        for (int z = 0; z < 9; z++)
          deltaF(z) = pFpu(z, y);
        
        MATRIX3 dF;
        MATRIX_UTIL::repackVec9(deltaF, dF);
        
        MATRIX3 dFhat = R.transpose() * dF;
        
        MATRIX3 deltaP = material->rotatedFirstPKDifferential(L, dFhat); 
              
        // rotate deltaP back
        deltaP = R * deltaP;

        VEC3F forceVecs[4];
        for(int z = 0; z < 4; z++){
          forceVecs[z] = deltaP * b[z];
        }
     
        // copy result into stiffness column
        for (int z = 0; z < 4; z++){
          if(indices[z] >= _tetMesh->unconstrainedNodes())
            continue;
          int row = indices[z] * 3;
          for (int a = 0; a < 3; a++){
            tetStiffness(z * 3 + a, y) = -forceVecs[z][a];
          }
        }
      }
      if(_tetMesh->isSkinningBasis()){
        _KCopies[id] += keyTetWeights[x] * _keyTetUs[x].transpose() * _keyTetTransforms[x].transpose() * tetStiffness * _keyTetTransforms[x] * _keyTetUs[x];
      }else{
        _KCopies[id] += keyTetWeights[x] * _keyTetUs[x].transpose() * tetStiffness * _keyTetUs[x];
      }
    }

  } // OMP

  for(int x = 0; x < _tetMesh->totalCores(); x++)
    stiffness += _KCopies[x];

  TIMING_BREAKDOWN::toc("Compute Reduced Stiffness Matrix");

  return stiffness;
}
