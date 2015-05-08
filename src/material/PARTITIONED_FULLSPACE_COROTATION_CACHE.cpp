#include <material/PARTITIONED_FULLSPACE_COROTATION_CACHE.h>
#include <util/MATRIX_UTIL.h>
#include <util/TIMING_BREAKDOWN.h>
#include <set>

#if USING_OPENMP
#include <omp.h>
#endif

PARTITIONED_FULLSPACE_COROTATION_CACHE::PARTITIONED_FULLSPACE_COROTATION_CACHE(TET_MESH* tetMesh):
  _tetMesh(tetMesh)
{
  int totalTets = tetMesh->tets().size();
  _internalForces.resize(_tetMesh->totalPartitions());
  _Rs.resize(totalTets);
  _Ss.resize(totalTets);
  _Ls.resize(totalTets);
}
PARTITIONED_FULLSPACE_COROTATION_CACHE::~PARTITIONED_FULLSPACE_COROTATION_CACHE()
{
}

void PARTITIONED_FULLSPACE_COROTATION_CACHE::cacheDecompositions()
{
  TIMING_BREAKDOWN::tic();

  vector<TET>& tets = _tetMesh->tets();
  _tetMesh->updateFullMesh();
  _tetMesh->generateF();
  vector<MATRIX3>& Fs = _tetMesh->F();

#if USING_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_OPENMP
    const int id = omp_get_thread_num();
    #else
    const int id = 0;
    #endif

    #if USING_OPENMP
    #pragma omp for schedule(static)
    #endif
    for(int x = 0; x < tets.size(); x++){
      int materialIndex = tets[x].materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
      MATRIX3& F = Fs[x];
      MATRIX3 U, V, Fhat;
      material->decomposeF(F, U, Fhat, V, _Rs[x], _Ss[x], _Ls[x]);
    }
  }
  TIMING_BREAKDOWN::toc("Cache Material Decompositions");
}

void PARTITIONED_FULLSPACE_COROTATION_CACHE::cachePartialDecompositions()
{
  TIMING_BREAKDOWN::tic();

  vector<TET>& tets = _tetMesh->tets();
  _tetMesh->updateFullMesh();

#if USING_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_OPENMP
    const int id = omp_get_thread_num();
    #else
    const int id = 0;
    #endif

      
    for(int x = 0; x < _tetMesh->totalPartitions(); x++){
      if(_tetMesh->partitionedFullsimDofs(x) == 0)
        continue;
      vector<int>& tetIDs = _tetMesh->partitionedFullsimTetIDs(x);
      // vector<int>& tetIDs = _tetMesh->partitionedTets(x);
      #if USING_OPENMP
      #pragma omp for schedule(static)
      #endif
      for(unsigned int y = 0; y < tetIDs.size(); y++){
        int tetID = tetIDs[y];

        int materialIndex = tets[tetID].materialIndex();
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
        MATRIX3 F = tets[tetID].F();
        MATRIX3 U, V, Fhat;
        material->decomposeF(F, U, Fhat, V, _Rs[tetID], _Ss[tetID], _Ls[tetID]);
      }
    }
    
  }
  TIMING_BREAKDOWN::toc("Cache Material Decompositions");
}


VECTOR& PARTITIONED_FULLSPACE_COROTATION_CACHE::computeInternalForce()
{
  TIMING_BREAKDOWN::tic();

  VECTOR& R = _tetMesh->partitionedInternalForce();
  if(R.size() != _tetMesh->partitionedDofs())
    R.resize(_tetMesh->partitionedDofs());

  R.setZero();

  vector<TET>& tets = _tetMesh->tets();

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif

    #if USING_OPENMP
    #pragma omp for schedule(static)
    #endif
    for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
      vector<int>& tetIDs = _tetMesh->partitionedTets(x);
      map<int, int>& partitionedVertexIDs = _tetMesh->partitionedVertexIDs(x);
      for(unsigned int y = 0; y < tetIDs.size(); y++){
        int tetID = tetIDs[y];
        TET& tet = tets[tetID];
        // compute the forces
        VECTOR forces(12);
        const VEC3F* b = tet.b();

        int materialIndex = tet.materialIndex();
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

        MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[tetID], _Ss[tetID]);
        forces.segment<3>(0) = firstPK * b[0];
        forces.segment<3>(3) = firstPK * b[1];
        forces.segment<3>(6) = firstPK * b[2];
        forces.segment<3>(9) = firstPK * b[3];

        for(int z = 0; z < 4; z++){
          int originalVertexID = _tetMesh->vertexID(tet.vertices[z]);
          if(!_tetMesh->isConstrained(originalVertexID)){
            int index = partitionedVertexIDs[originalVertexID];
            R.segment<3>(_tetMesh->partitionDofStartIdx(x) + index * 3) += forces.segment<3>(z * 3);
          }
        }
      }
    }
  }

  TIMING_BREAKDOWN::toc("Compute Partitioned Internal Force");
  return R;
}

Real PARTITIONED_FULLSPACE_COROTATION_CACHE::computeElasticEnergy()
{
  TIMING_BREAKDOWN::tic();

  Real elasticEnergy = 0.0;

  vector<TET>& tets = _tetMesh->tets();
#if USING_OPENMP
#pragma omp parallel
#endif
  {
    #if USING_OPENMP
    const int id = omp_get_thread_num();
    #else
    const int id = 0;
    #endif

    #if USING_OPENMP
    #pragma omp for schedule(static) reduction(+ : elasticEnergy)
    #endif
    for (int x = 0; x < tets.size(); x++)
    {
      int materialIndex = tets[x].materialIndex();
      
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      elasticEnergy += tets[x].restVolume() * material->strainEnergy(_Rs[x], _Ss[x]);
    }
  }
  TIMING_BREAKDOWN::toc("Compute Elastic Energy");
  return elasticEnergy;
}

BLOCK_COO_MATRIX& PARTITIONED_FULLSPACE_COROTATION_CACHE::computeStiffnessMatrix()
{
  // TIMING_BREAKDOWN::tic();

  BLOCK_COO_MATRIX& stiffness = _tetMesh->partitionedStiffnessMatrix();
  stiffness.resizeAndWipe(_tetMesh->totalPartitions(), _tetMesh->totalPartitions());

  vector<int> rowSizes(_tetMesh->totalPartitions());
  vector<int> colSizes(_tetMesh->totalPartitions());
  
  for(int x = 0; x < _tetMesh->totalPartitions(); x++){
    rowSizes[x] = colSizes[x] = _tetMesh->partitionedVertices(x).size() * 3;
  }
  stiffness.setBlockDimensions(rowSizes, colSizes);

  vector<TET>& tets = _tetMesh->tets();
  

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif

    #if USING_OPENMP
    #pragma omp for schedule(static)
    #endif
    for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
      
      vector<int>& tetIDs = _tetMesh->partitionedTets(x);
      map<int, int>& partitionedVertexIDs = _tetMesh->partitionedVertexIDs(x);

      int dofs = _tetMesh->partitionedVertices(x).size() * 3;
      COO_MATRIX& localCooStiffness = stiffness(x, x);

      // cout << " partition " << x << " stiffness dimension " << localCooStiffness.rows() << " " << localCooStiffness.cols() << endl;

      for(unsigned int t = 0; t < tetIDs.size(); t++){
        int tetID = tetIDs[t];
        TET& tet = tets[tetID];

        int materialIndex = tet.materialIndex();

        const VEC3F* b = tet.b();

        int indices[] = {_tetMesh->vertexID(tet.vertices[0]),
                       _tetMesh->vertexID(tet.vertices[1]),
                       _tetMesh->vertexID(tet.vertices[2]),
                       _tetMesh->vertexID(tet.vertices[3])};

        MATRIX& pFpu = tet.PFPu();

        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
        
        MATRIX3& R = _Rs[tetID];
        MATRIX3& S = _Ss[tetID];
        MATRIX3& L = _Ls[tetID];
        // compute each column of the matrix
        for (int y = 0; y < 12; y++)
        {
          int col = indices[y / 3] * 3 + y % 3;
          if(col >= _tetMesh->dofs())
            continue;

          col = partitionedVertexIDs[indices[y / 3]] * 3 + y % 3;
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
       
          // copy result into stiffness column
          for (int z = 0; z < 4; z++){
            if(indices[z] >= _tetMesh->unconstrainedNodes())
              continue;
            int row = partitionedVertexIDs[indices[z]] * 3;
            for (int a = 0; a < 3; a++){
              localCooStiffness.add(-forceVecs(z * 3 + a), row + a, col);
            }
          }
        }
      }
      localCooStiffness.order();
      localCooStiffness.aggregate();
    }
  }

  // TIMING_BREAKDOWN::toc("Compute Partitioned Stiffness Matrix");

  return stiffness;
}

VECTOR& PARTITIONED_FULLSPACE_COROTATION_CACHE::computePartialInternalForce(int partition)
{
  int dofs = _tetMesh->partitionedVertices(partition).size() * 3;
  _internalForces[partition].conservativeResize(dofs);
  _internalForces[partition].setZero();

  vector<TET>& tets = _tetMesh->tets();
  // vector<int>& tetIDs = _tetMesh->partitionedTets(partition);
  vector<int>& tetIDs = _tetMesh->partitionedFullsimTetIDs(partition);
  
  map<int, int>& partitionedVertexIDs = _tetMesh->partitionedVertexIDs(partition);
  vector<int>& newOrderings = _tetMesh->adaptivePartitionedVertexOrdering(partition);

#if USING_OPENMP
#pragma omp parallel
#endif
  {
    VECTOR localR(dofs);
    localR.setZero();

#if USING_OPENMP
    const int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
    const int id  = 0;
#endif
    for(unsigned int x = 0; x < tetIDs.size(); x++){
      int tetID = tetIDs[x];
      TET& tet = tets[tetID];
      // compute the forces
      // VECTOR forces(12);
      // const VEC3F* b = tet.b();

      int materialIndex = tet.materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[tetID], _Ss[tetID]);

      MATRIX3x4 forces = firstPK * tet.bMat();
      // forces.segment<3>(0) = firstPK * b[0];
      // forces.segment<3>(3) = firstPK * b[1];
      // forces.segment<3>(6) = firstPK * b[2];
      // forces.segment<3>(9) = firstPK * b[3];

      for(int z = 0; z < 4; z++){
        int originalVertexID = _tetMesh->vertexID(tet.vertices[z]);
        if(_tetMesh->isConstrained(originalVertexID))
          continue;

        int index = newOrderings[partitionedVertexIDs[originalVertexID]];

        localR.segment<3>(index * 3) += forces.col(z);//forces.segment<3>(z * 3);
      }
    }
    #if USING_OPENMP
    #pragma omp critical
    #endif
    {
       _internalForces[partition] -= localR;
    }
  }
  return _internalForces[partition];
}
vector<VECTOR>& PARTITIONED_FULLSPACE_COROTATION_CACHE::computePartialInternalForce()
{
  TIMING_BREAKDOWN::tic();

  for(unsigned int x = 0; x < _tetMesh->totalPartitions(); x++){
    if(_tetMesh->partitionedFullsimDofs(x) > 0)
      computePartialInternalForce(x);
  }
  
  TIMING_BREAKDOWN::toc("Compute Partial Partitioned Internal Force");
  return _internalForces;
}
void PARTITIONED_FULLSPACE_COROTATION_CACHE::computePartialStiffnessMatrix(int partition, COO_MATRIX& diagMat, COO_MATRIX& offDiagMat, COO_MATRIX& diagMat2)
{
  vector<TET>& tets = _tetMesh->tets();
  int fullDofs = _tetMesh->partitionedFullsimDofs(partition);
  int reducedDofs = _tetMesh->partitionedReducedsimDofs(partition);

  diagMat.resize(fullDofs, fullDofs);
  offDiagMat.resize(reducedDofs, fullDofs);
  diagMat2.resize(reducedDofs, reducedDofs);

  vector<int>& fullsimTetIDs = _tetMesh->partitionedFullsimTetIDs(partition);
  map<int, int>& partitionedVertexIDs = _tetMesh->partitionedVertexIDs(partition);

  vector<int>& newOrderings = _tetMesh->adaptivePartitionedVertexOrdering(partition);

  vector<int> diagOffsets(_tetMesh->totalCores() + 1, 0);
  vector<int> offDiagOffsets(_tetMesh->totalCores() + 1, 0);
  vector<int> diag2Offsets(_tetMesh->totalCores() + 1, 0);

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif

  COO_MATRIX diagLocalCooStiffness(fullDofs, fullDofs);
  COO_MATRIX offDiagLocalCooStiffness(reducedDofs, fullDofs);
  COO_MATRIX diag2LocalCooStiffness(reducedDofs, reducedDofs);

  #if USING_OPENMP
  #pragma omp for schedule(static)
  #endif
  for (unsigned int x = 0; x < fullsimTetIDs.size(); x++)
  {
    int tetID = fullsimTetIDs[x];
    TET& tet = tets[tetID];
    int materialIndex = tet.materialIndex();
    // compute Bm matrix
    // const VEC3F* b = tet.b();

    int indices[] = {_tetMesh->vertexID(tet.vertices[0]),
                     _tetMesh->vertexID(tet.vertices[1]),
                     _tetMesh->vertexID(tet.vertices[2]),
                     _tetMesh->vertexID(tet.vertices[3])};
    
    // get PFPu from the tet - only depends on rest state,
    // so don't need to update tet state at all
    MATRIX& pFpu = tet.PFPu();
    
    COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
    
    MATRIX3& R = _Rs[tetID];
    MATRIX3& S = _Ss[tetID];
    MATRIX3& L = _Ls[tetID];
    // compute each column of the matrix
    for (int y = 0; y < 12; y++)
    {
      if(_tetMesh->isConstrained(indices[y / 3]))
        continue;

      int partitionedVertexID = partitionedVertexIDs[indices[y / 3]];
      partitionedVertexID = newOrderings[partitionedVertexID];

      // if(partitionedVertexID >= fullDofs / 3)
          // continue;
      int col = partitionedVertexID * 3 + y % 3;

      // extract a column from pFpu (ie pretend 'x' is just a 
      // Kronecker delta)
      VECTOR deltaF = pFpu.col(y);
      
      MATRIX3 dF;
      MATRIX_UTIL::repackVec9(deltaF, dF);
      
      MATRIX3 dFhat = R.transpose() * dF;
      
      MATRIX3 deltaP = material->rotatedFirstPKDifferential(L, dFhat); 
            
      // rotate deltaP back
      deltaP = R * deltaP;

      // VECTOR forceVecs(12);
      // for(int z = 0; z < 4; z++){
        // forceVecs.segment<3>(z * 3) = deltaP * b[z];
      // }
      MATRIX3x4 forceVecs = deltaP * tet.bMat();
   
      if(col < fullDofs){
        // copy result into stiffness column
        for (int z = 0; z < 4; z++){
          if(_tetMesh->isConstrained(indices[z]))
            continue;
          
          int vid = partitionedVertexIDs[indices[z]];
          vid = newOrderings[vid];

          int row = vid * 3;
          if(row < fullDofs){
            for(int a = 0; a < 3; a++)
              diagLocalCooStiffness.add(-forceVecs(a, z), row + a, col);
          }else{
            row -= fullDofs;
            for(int a = 0; a < 3; a++)
              offDiagLocalCooStiffness.add(-forceVecs(a, z), row + a, col);
          }
        }
      }else{
        col -= fullDofs;
        // copy result into stiffness column
        for (int z = 0; z < 4; z++){
          if(_tetMesh->isConstrained(indices[z]))
            continue;
          
          int vid = partitionedVertexIDs[indices[z]];
          vid = newOrderings[vid];

          int row = vid * 3;
          if(row >= fullDofs){
            row -= fullDofs;
            for(int a = 0; a < 3; a++)
              diag2LocalCooStiffness.add(-forceVecs(a, z), row + a, col);
          }
        }
      }
    }
       
  }

  diagLocalCooStiffness.order();
  diagLocalCooStiffness.aggregate();

  offDiagLocalCooStiffness.order();
  offDiagLocalCooStiffness.aggregate();

  diag2LocalCooStiffness.order();
  diag2LocalCooStiffness.aggregate();

  diagOffsets[id + 1] = diagLocalCooStiffness.nnZ();
  offDiagOffsets[id + 1] = offDiagLocalCooStiffness.nnZ();
  diag2Offsets[id + 1] = diag2LocalCooStiffness.nnZ();

#if USING_OPENMP
#pragma omp barrier // ===========================
#pragma omp single
#endif
  { 
      for (int t = 1; t < diagOffsets.size(); t++){
        diagOffsets[t] += diagOffsets[t-1];
        offDiagOffsets[t] += offDiagOffsets[t - 1];
        diag2Offsets[t] += diag2Offsets[t-1];
      }
      diagMat.matrix().resize(diagOffsets.back());
      offDiagMat.matrix().resize(offDiagOffsets.back());
      diagMat2.matrix().resize(diag2Offsets.back());
  }
  vector<TRIPLET>::iterator diagDest = diagMat.matrix().begin();
  std::advance(diagDest, diagOffsets[id]);
  std::copy(diagLocalCooStiffness.matrix().begin(), diagLocalCooStiffness.matrix().end(), diagDest);

  vector<TRIPLET>::iterator offDiagDest = offDiagMat.matrix().begin();
  std::advance(offDiagDest, offDiagOffsets[id]);
  std::copy(offDiagLocalCooStiffness.matrix().begin(), offDiagLocalCooStiffness.matrix().end(), offDiagDest);

  vector<TRIPLET>::iterator diagDest2 = diagMat2.matrix().begin();
  std::advance(diagDest2, diag2Offsets[id]);
  std::copy(diag2LocalCooStiffness.matrix().begin(), diag2LocalCooStiffness.matrix().end(), diagDest2);

  } // OMP
}
