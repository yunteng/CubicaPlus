#include <material/PARTITIONED_SUBSPACE_COROTATION_CACHE.h>
#include <util/MATRIX_UTIL.h>
#include <util/TIMING_BREAKDOWN.h>
#include <set>

#if USING_SUBSPACE_OPENMP
#include <omp.h>
#endif

PARTITIONED_SUBSPACE_COROTATION_CACHE::PARTITIONED_SUBSPACE_COROTATION_CACHE(SUBSPACE_TET_MESH* tetMesh):
  _tetMesh(tetMesh),
  _coupledWithFullspace(false)
{
  int totalKeyTets = _tetMesh->keyTets().size();
  _Rs.resize(totalKeyTets);
  _Ss.resize(totalKeyTets);
  _Ls.resize(totalKeyTets);
  _keyTetTransforms.resize(totalKeyTets);

  for(unsigned int x = 0; x < totalKeyTets; x++){
    _keyTetTransforms[x].resize(12, 12);
    _keyTetTransforms[x].reserve(Eigen::VectorXi::Constant(12, 3));
  }

  _keyTetUs.resize(tetMesh->totalPartitions());

  vector<TET>& tets = _tetMesh->tets();

  // compute the keyTetUs
  for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
    vector<int>& keyTetsIDs = _tetMesh->partitionedKeyTets(x);
    for(unsigned int y = 0; y < keyTetsIDs.size(); y++){
      int tetID = keyTetsIDs[y];
      TET& tet = tets[tetID];
      MATRIX U(12, _tetMesh->partitionRank(x));
      U.setZero();

      for(int z = 0; z < 4; z++){
        int vertexID = _tetMesh->vertexID(tet.vertices[z]);
        if(_tetMesh->isConstrained(vertexID))
          continue;
        int partitionedVertexID = _tetMesh->partitionedVertexID(x, vertexID);
        U.block(z * 3, 0, 3, _tetMesh->partitionRank(x)) = _tetMesh->partitionVertexBasis(x, partitionedVertexID);
      }
      _keyTetUs[x].push_back(U);
    }
  }
  _isFullsimTets.resize(_tetMesh->tets().size());
  _isFullsimTets.setZero();
}
PARTITIONED_SUBSPACE_COROTATION_CACHE::~PARTITIONED_SUBSPACE_COROTATION_CACHE()
{
}

void PARTITIONED_SUBSPACE_COROTATION_CACHE::cacheDecompositions()
{
  TIMING_BREAKDOWN::tic();

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& keyTetsIDs = _tetMesh->keyTets();
  /* already done by the fullspace corotation cache*/
  if(!_coupledWithFullspace)
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

void PARTITIONED_SUBSPACE_COROTATION_CACHE::cacheKeyTetTransforms(vector<MATRIX3>& transformMatrix)
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
void PARTITIONED_SUBSPACE_COROTATION_CACHE::registerFullsimTets()
{
  _isFullsimTets.setZero();
  #if USING_SUBSPACE_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_tetMesh->partitionedFullsimDofs(x) == 0)
      continue;
    vector<int>& tetIDs = _tetMesh->partitionedFullsimTetIDs(x);
    for(unsigned int y = 0; y < tetIDs.size(); y++)
      _isFullsimTets[tetIDs[y]] = 1;
  }
}

VECTOR& PARTITIONED_SUBSPACE_COROTATION_CACHE::computeInternalForce()
{
  TIMING_BREAKDOWN::tic();

  VECTOR& R = _tetMesh->reducedInternalForce();
  if(R.size() != _tetMesh->totalPartitionRank())
    R.resize(_tetMesh->totalPartitionRank());

  R.setZero();

  vector<TET>& tets = _tetMesh->tets();

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_SUBSPACE_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif

    #if USING_SUBSPACE_OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
      // if(_coupledWithFullspace && _tetMesh->partitionedFullsimDofs(x) > 0)
        // continue;

      vector<int>& keyTetsIDs = _tetMesh->partitionedKeyTets(x);
      vector<Real>& keyTetWeights = _tetMesh->partitionedKeyWeights(x);
      for(unsigned int y = 0; y < keyTetsIDs.size(); y++){
        int index = _tetMesh->partitionKeyTetStartIdx(x) + y;
        int tetID = keyTetsIDs[y];
        if(_isFullsimTets[tetID] > 0)
          continue;
        // compute the forces
        VECTOR forces(12);
        const VEC3F* b = tets[tetID].b();

        int materialIndex = tets[tetID].materialIndex();
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

        MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[index], _Ss[index]);
        forces.segment<3>(0) = firstPK * b[0];
        forces.segment<3>(3) = firstPK * b[1];
        forces.segment<3>(6) = firstPK * b[2];
        forces.segment<3>(9) = firstPK * b[3];
        R.segment(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x)) +=  _keyTetUs[x][y].transpose() * (_keyTetTransforms[index].transpose() * (keyTetWeights[y] * forces));
      }
    }
  }
  TIMING_BREAKDOWN::toc("Compute Reduced Internal Force");
  return R;
}

Real PARTITIONED_SUBSPACE_COROTATION_CACHE::computeElasticEnergy()
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

MATRIX& PARTITIONED_SUBSPACE_COROTATION_CACHE::computeReducedStiffnessMatrix()
{
  TIMING_BREAKDOWN::tic();

  MATRIX& stiffness = _tetMesh->reducedStiffnessMatrix();

  if(stiffness.rows() != _tetMesh->totalPartitionRank())
    stiffness.resize(_tetMesh->totalPartitionRank(), _tetMesh->totalPartitionRank());
  stiffness.setZero();

  vector<TET>& tets = _tetMesh->tets();
  

#if USING_SUBSPACE_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_SUBSPACE_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif

    #if USING_SUBSPACE_OPENMP
    #pragma omp for schedule(dynamic)
    #endif
    for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
      // if(_coupledWithFullspace && _tetMesh->partitionedFullsimDofs(x) > 0)
        // continue;

      vector<int>& keyTetsIDs = _tetMesh->partitionedKeyTets(x);
      vector<Real>& keyTetWeights = _tetMesh->partitionedKeyWeights(x);
      for(unsigned int t = 0; t < keyTetsIDs.size(); t++){
        MATRIX tetStiffness(12, 12);
        tetStiffness.setZero();

        int index = _tetMesh->partitionKeyTetStartIdx(x) + t;
        int tetID = keyTetsIDs[t];
        if(_isFullsimTets[tetID] > 0)
          continue;
        
        TET& tet = tets[tetID];

        const VEC3F* b = tets[tetID].b();

        int materialIndex = tets[tetID].materialIndex();
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

        int indices[] = {_tetMesh->vertexID(tet.vertices[0]),
                         _tetMesh->vertexID(tet.vertices[1]),
                         _tetMesh->vertexID(tet.vertices[2]),
                         _tetMesh->vertexID(tet.vertices[3])};

        MATRIX& pFpu = tet.PFPu();

        MATRIX3& R = _Rs[index];
        MATRIX3& S = _Ss[index];
        MATRIX3& L = _Ls[index];
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

          VECTOR forceVecs(12);
          for(int z = 0; z < 4; z++){
            forceVecs.segment<3>(z * 3) = deltaP * b[z];
          }
          // apply key tet weights!!!
          forceVecs *= keyTetWeights[t];
       
          // copy result into stiffness column
          for (int z = 0; z < 4; z++){
            if(indices[z] >= _tetMesh->unconstrainedNodes())
              continue;
            int row = indices[z] * 3;
            for (int a = 0; a < 3; a++){
              tetStiffness(z * 3 + a, y) = -forceVecs(z * 3 + a);
            }
          }
        }
        stiffness.block(_tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRankStartIdx(x), _tetMesh->partitionRank(x), _tetMesh->partitionRank(x)) += _keyTetUs[x][t].transpose() * _keyTetTransforms[index].transpose() * tetStiffness * _keyTetTransforms[index] * _keyTetUs[x][t];
      }
    }
  }

  TIMING_BREAKDOWN::toc("Compute Reduced Stiffness Matrix");

  return stiffness;
}

