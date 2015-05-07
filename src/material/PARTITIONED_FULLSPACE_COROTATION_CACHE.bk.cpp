void PARTITIONED_FULLSPACE_COROTATION_CACHE::computePartialStiffnessMatrix(int partition, COO_MATRIX& diagMat, COO_MATRIX& offDiagMat, COO_MATRIX& internalMat)
{
  vector<TET>& tets = _tetMesh->tets();
  int fullDofs = _tetMesh->partitionedFullsimDofs(partition);
  int reducedDofs = _tetMesh->partitionedReducedsimDofs(partition);

  diagMat.resize(fullDofs, fullDofs);
  offDiagMat.resize(reducedDofs, fullDofs);
  internalMat.resize(reducedDofs, reducedDofs);

  vector<int>& fullsimTetIDs = _tetMesh->partitionedTets(partition);
  map<int, int>& partitionedVertexIDs = _tetMesh->partitionedVertexIDs(partition);

  vector<int>& newOrderings = _tetMesh->adaptivePartitionedVertexOrdering(partition);

  vector<int> diagOffsets(_tetMesh->totalCores() + 1, 0);
  vector<int> offDiagOffsets(_tetMesh->totalCores() + 1, 0);
  vector<int> internalOffsets(_tetMesh->totalCores() + 1, 0);

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
  COO_MATRIX internalLocalCooStiffness(reducedDofs, reducedDofs);

  #if USING_OPENMP
  #pragma omp for schedule(static)
  #endif
  for (unsigned int x = 0; x < fullsimTetIDs.size(); x++)
  {
    int tetID = fullsimTetIDs[x];
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

      VECTOR forceVecs(12);
      for(int z = 0; z < 4; z++){
        forceVecs.segment<3>(z * 3) = deltaP * b[z];
      }
   
      // copy result into stiffness column
      for (int z = 0; z < 4; z++){
        if(_tetMesh->isConstrained(indices[z]))
          continue;
        
        int vid = partitionedVertexIDs[indices[z]];
        vid = newOrderings[vid];

        int row = vid * 3;
        if(col < fullDofs){
          if(row < fullDofs){
            for(int a = 0; a < 3; a++)
              diagLocalCooStiffness.add(-forceVecs(z * 3 + a), row + a, col);
          }else{
            row -= fullDofs;
            for(int a = 0; a < 3; a++)
              offDiagLocalCooStiffness.add(-forceVecs(z * 3 + a), row + a, col);
          }
        }else{
          if(row >= fullDofs){
            row -= fullDofs;
            for(int a = 0; a < 3; a++)
              internalLocalCooStiffness.add(-forceVecs(z * 3 + a), row + a, col - fullDofs);
          }
        }
      }
    }
       
  }

  diagLocalCooStiffness.order();
  diagLocalCooStiffness.aggregate();

  offDiagLocalCooStiffness.order();
  offDiagLocalCooStiffness.aggregate();

  internalLocalCooStiffness.order();
  internalLocalCooStiffness.aggregate();

  diagOffsets[id + 1] = diagLocalCooStiffness.nnZ();
  offDiagOffsets[id + 1] = offDiagLocalCooStiffness.nnZ();
  internalOffsets[id + 1] = internalLocalCooStiffness.nnZ();

#if USING_OPENMP
#pragma omp barrier // ===========================
#pragma omp single
#endif
  { 
      for (int t = 1; t < diagOffsets.size(); t++){
        diagOffsets[t] += diagOffsets[t-1];
        offDiagOffsets[t] += offDiagOffsets[t - 1];
        internalOffsets[t] += internalOffsets[t - 1];
      }
      diagMat.matrix().resize(diagOffsets.back());
      offDiagMat.matrix().resize(offDiagOffsets.back());
      internalMat.matrix().resize(internalOffsets.back());
  }
  vector<TRIPLET>::iterator diagDest = diagMat.matrix().begin();
  std::advance(diagDest, diagOffsets[id]);
  std::copy(diagLocalCooStiffness.matrix().begin(), diagLocalCooStiffness.matrix().end(), diagDest);

  vector<TRIPLET>::iterator offDiagDest = offDiagMat.matrix().begin();
  std::advance(offDiagDest, offDiagOffsets[id]);
  std::copy(offDiagLocalCooStiffness.matrix().begin(), offDiagLocalCooStiffness.matrix().end(), offDiagDest);

  vector<TRIPLET>::iterator internalDest = internalMat.matrix().begin();
  std::advance(internalDest, internalOffsets[id]);
  std::copy(internalLocalCooStiffness.matrix().begin(), internalLocalCooStiffness.matrix().end(), internalDest);
  } // OMP
}
