#include <material/FULLSPACE_COROTATION_CACHE.h>
#include <util/MATRIX_UTIL.h>
#include <util/TIMING_BREAKDOWN.h>

FULLSPACE_COROTATION_CACHE::FULLSPACE_COROTATION_CACHE(TET_MESH* tetMesh):
  _tetMesh(tetMesh)
{
  int totalTets = tetMesh->tets().size();
  _Rs.resize(totalTets);
  _Ss.resize(totalTets);
  _Ls.resize(totalTets);

  int rank = tetMesh->dofs();
  _RCopies = new VECTOR[tetMesh->totalCores()];
  for(int x = 0; x < tetMesh->totalCores(); x++)
    _RCopies[x].resize(rank);
}
FULLSPACE_COROTATION_CACHE::~FULLSPACE_COROTATION_CACHE()
{
  delete[] _RCopies;
}

void FULLSPACE_COROTATION_CACHE::cacheDecompositions(bool surfaceOnly)
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

    if(!surfaceOnly){
      #if USING_OPENMP
      #pragma omp for schedule(static)
      #endif
      for(int x = 0; x < tets.size(); x++){
        int materialIndex = tets[x].materialIndex();
        // cout << x << " material index " << materialIndex << endl;
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
        MATRIX3& F = Fs[x];
        MATRIX3 U, V, Fhat;
        material->decomposeF(F, U, Fhat, V, _Rs[x], _Ss[x], _Ls[x]);
      }
    }else{
      vector<int>& fullsimTetIDs = _tetMesh->surfaceTetIDs();
      #if USING_OPENMP
      #pragma omp for schedule(static)
      #endif
      for(int x = 0; x < fullsimTetIDs.size(); x++){
        int tetID = fullsimTetIDs[x];
        int materialIndex = tets[tetID].materialIndex();
        COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);
        MATRIX3& F = Fs[tetID];
        MATRIX3 U, V, Fhat;
        material->decomposeF(F, U, Fhat, V, _Rs[tetID], _Ss[tetID], _Ls[tetID]);
      }
    }
  }
  TIMING_BREAKDOWN::toc("Cache Material Decompositions"); 
}
void FULLSPACE_COROTATION_CACHE::computeIndividualTetInternalForces(VECTOR& output)
{
  vector<TET>& tets = _tetMesh->tets();
  output.resize(tets.size() * 12);

  for(unsigned int x = 0; x < tets.size(); x++)
  {
    // compute the forces
    VECTOR forces(12);
    const VEC3F* b = tets[x].b();

    int materialIndex = tets[x].materialIndex();
    COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[0][materialIndex]);

    MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[x], _Ss[x]);
    // forces[0] = firstPK * b[0];
    // forces[1] = firstPK * b[1];
    // forces[2] = firstPK * b[2];
    // forces[3] = firstPK * b[3];

    forces.segment<3>(0) = firstPK * b[0];
    forces.segment<3>(3) = firstPK * b[1];
    forces.segment<3>(6) = firstPK * b[2];
    forces.segment<3>(9) = firstPK * b[3];

    output.segment<12>(x * 12) = forces;
  }
}

void FULLSPACE_COROTATION_CACHE::computeIndividualTetInternalForces(vector<int>& tetIDs, VECTOR& output)
{
  vector<TET>& tets = _tetMesh->tets();
  output.resize(tetIDs.size() * 12);

  for(unsigned int x = 0; x < tetIDs.size(); x++)
  {
    int id = tetIDs[x];

    VECTOR forces(12);

    const VEC3F* b = tets[id].b();

    int materialIndex = tets[id].materialIndex();
    COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[0][materialIndex]);

    MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[id], _Ss[id]);
    // forces[0] = firstPK * b[0];
    // forces[1] = firstPK * b[1];
    // forces[2] = firstPK * b[2];
    // forces[3] = firstPK * b[3];

    forces.segment<3>(0) = firstPK * b[0];
    forces.segment<3>(3) = firstPK * b[1];
    forces.segment<3>(6) = firstPK * b[2];
    forces.segment<3>(9) = firstPK * b[3];

    output.segment<12>(x * 12) = forces;
  }
}

VECTOR& FULLSPACE_COROTATION_CACHE::computeInternalForce()
{
  TIMING_BREAKDOWN::tic();

  VECTOR& R = _tetMesh->internalForce();
  if(R.size() != _tetMesh->dofs())
    R.resize(_tetMesh->dofs());
  
  R.setZero();

  for (int x = 0; x < _tetMesh->totalCores(); x++){
    _RCopies[x].setZero();
  }

  vector<TET>& tets = _tetMesh->tets();

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
    const int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
    const int id  = 0;
#endif
  // populate the forces
    for (unsigned int x = 0; x < tets.size(); x++)
    {
      // compute the forces
      VEC3F forces[4];
      const VEC3F* b = tets[x].b();

      int materialIndex = tets[x].materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[x], _Ss[x]);
      forces[0] = firstPK * b[0];
      forces[1] = firstPK * b[1];
      forces[2] = firstPK * b[2];
      forces[3] = firstPK * b[3];

      // add forces to corresponding vertices
      for (int y = 0; y < 4; y++)
      {
        int index = _tetMesh->vertexID(tets[x].vertices[y]);
        if (!_tetMesh->isConstrained(index))
        {
          _RCopies[id].segment<3>(index * 3) += forces[y];
        }
      }
    }
  }
  // merge all the copies
  for (int x = 0; x < _tetMesh->totalCores(); x++)
    R += _RCopies[x];

  TIMING_BREAKDOWN::toc("Compute Internal Force");
  return R;
}

Real FULLSPACE_COROTATION_CACHE::computeElasticEnergy()
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

COO_MATRIX& FULLSPACE_COROTATION_CACHE::computeStiffnessMatrix()
{
  TIMING_BREAKDOWN::tic();

  COO_MATRIX& stiffness = _tetMesh->stiffnessMatrix();

  if(stiffness.rows() != _tetMesh->dofs())
    stiffness.resize(_tetMesh->dofs(), _tetMesh->dofs());

  vector<TET>& tets = _tetMesh->tets();

  vector<int> offsets(_tetMesh->totalCores() + 1, 0);

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
    #if USING_OPENMP
    const int id  = omp_get_thread_num();
    #else
    const int id = 0;
    #endif
    
    COO_MATRIX localCooStiffness(_tetMesh->dofs(), _tetMesh->dofs());

    // vector<pair<pair<int, int>, Real> >& mat = localCooStiffness.matrix();

    #if USING_OPENMP
    #pragma omp for  schedule(static)
    #endif
    for (unsigned int x = 0; x < tets.size(); x++)
    {
      TET& tet = tets[x];
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
            // if(row + a >= col)
              // mat.push_back(make_pair(make_pair(row + a, col), -forceVecs[z][a]));
            localCooStiffness.add(-forceVecs[z][a], row + a, col);
          }
        }
      }
         
    }

    localCooStiffness.order();
    localCooStiffness.aggregate();
    offsets[id + 1] = localCooStiffness.nnZ();

#if USING_OPENMP
#pragma omp barrier // ===========================
#pragma omp single
#endif
    { 
        for (int t = 1; t < offsets.size(); t++)
            offsets[t] += offsets[t-1];
        stiffness.matrix().resize(offsets.back());
    }

    vector<TRIPLET>::iterator dest = stiffness.matrix().begin();
    std::advance(dest, offsets[id]);
    std::copy(localCooStiffness.matrix().begin(), localCooStiffness.matrix().end(), dest);

  } // OMP
  stiffness.order();
  stiffness.aggregate(false);

  TIMING_BREAKDOWN::toc("Compute Stiffness Matrix");

  return stiffness;
}

BLOCK_SPARSE_MATRIX& FULLSPACE_COROTATION_CACHE::computeSurfaceStiffnessMatrix()
{
  TIMING_BREAKDOWN::tic();

  BLOCK_SPARSE_MATRIX& stiffness = _tetMesh->surfaceStiffnessMatrix();
  int surfaceSize = _tetMesh->surfaceVertexSize() * 3;
  int internalSize = _tetMesh->dofs() - surfaceSize;

  stiffness.resizeAndWipe(2, 1);

  vector<int> rowSizes(2);
  vector<int> colSizes(2);
  rowSizes[0] = colSizes[0] = surfaceSize;
  rowSizes[1] = internalSize;
  colSizes[1] = surfaceSize;

  stiffness.setBlockDimensions(rowSizes, colSizes);

  COO_MATRIX diagPart(surfaceSize, surfaceSize);

  COO_MATRIX offDiagPart(internalSize, surfaceSize);

  // if(stiffness.rows() != _tetMesh->dofs() || stiffness.cols() != surfaceSize)
    // stiffness.resize(_tetMesh->dofs(), surfaceSize);

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& surfaceTetIDs = _tetMesh->surfaceTetIDs();

  vector<int> diagOffsets(_tetMesh->totalCores() + 1, 0);
  vector<int> offDiagOffsets(_tetMesh->totalCores() + 1, 0);

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
    #if USING_OPENMP
    const int id  = omp_get_thread_num();
    #else
    const int id = 0;
    #endif
    
    COO_MATRIX diagLocalCooStiffness(surfaceSize, surfaceSize);
    COO_MATRIX offDiagLocalCooStiffness(internalSize, surfaceSize);

    // vector<pair<pair<int, int>, Real> >& mat = localCooStiffness.matrix();

    #if USING_OPENMP
    #pragma omp for  schedule(static)
    #endif
    for (unsigned int x = 0; x < surfaceTetIDs.size(); x++)
    {
      int tetID = surfaceTetIDs[x];
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
        int col = indices[y / 3] * 3 + y % 3;
        if(col >= surfaceSize)
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
          if(row < surfaceSize){
            for(int a = 0; a < 3; a++)
              diagLocalCooStiffness.add(-forceVecs[z][a], row + a, col);
          }else{
            row -= surfaceSize;
            for(int a = 0; a < 3; a++)
              offDiagLocalCooStiffness.add(-forceVecs[z][a], row + a, col);
          }
        }
      }
         
    }

    diagLocalCooStiffness.order();
    diagLocalCooStiffness.aggregate();

    offDiagLocalCooStiffness.order();
    offDiagLocalCooStiffness.aggregate();

    diagOffsets[id + 1] = diagLocalCooStiffness.nnZ();
    offDiagOffsets[id + 1] = offDiagLocalCooStiffness.nnZ();

#if USING_OPENMP
#pragma omp barrier // ===========================
#pragma omp single
#endif
    { 
        for (int t = 1; t < diagOffsets.size(); t++){
          diagOffsets[t] += diagOffsets[t-1];
          offDiagOffsets[t] += offDiagOffsets[t - 1];
        }
        diagPart.matrix().resize(diagOffsets.back());
        offDiagPart.matrix().resize(offDiagOffsets.back());
    }

    vector<TRIPLET>::iterator diagDest = diagPart.matrix().begin();
    std::advance(diagDest, diagOffsets[id]);
    std::copy(diagLocalCooStiffness.matrix().begin(), diagLocalCooStiffness.matrix().end(), diagDest);

    vector<TRIPLET>::iterator offDiagDest = offDiagPart.matrix().begin();
    std::advance(offDiagDest, offDiagOffsets[id]);
    std::copy(offDiagLocalCooStiffness.matrix().begin(), offDiagLocalCooStiffness.matrix().end(), offDiagDest);

  } // OMP
  // diagPart.order();
  // diagPart.aggregate(false);

  // offDiagPart.order();
  // offDiagPart.aggregate(false);
  stiffness.equals(diagPart, 0, 0);
  stiffness.equals(offDiagPart, 1, 0);

  TIMING_BREAKDOWN::toc("Compute Surface Stiffness Matrix");

  return stiffness;
}

VECTOR& FULLSPACE_COROTATION_CACHE::computeSurfaceInternalForce()
{
  TIMING_BREAKDOWN::tic();

  int surfaceSize = _tetMesh->surfaceVertexSize() * 3;
  VECTOR& R = _tetMesh->surfaceInternalForce();
  if(R.size() != surfaceSize)
    R.resize(surfaceSize);
  
  R.setZero();

  for (int x = 0; x < _tetMesh->totalCores(); x++){
    _RCopies[x].conservativeResize(surfaceSize);
    _RCopies[x].setZero();
  }

  vector<TET>& tets = _tetMesh->tets();
  vector<int>& surfaceTetIDs = _tetMesh->surfaceTetIDs();

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
#if USING_OPENMP
    const int id  = omp_get_thread_num();
#pragma omp for  schedule(static)
#else
    const int id  = 0;
#endif
  // populate the forces
    for (unsigned int x = 0; x < surfaceTetIDs.size(); x++)
    {
      int tetID = surfaceTetIDs[x];
      // compute the forces
      VEC3F forces[4];
      const VEC3F* b = tets[tetID].b();

      int materialIndex = tets[tetID].materialIndex();
      COROTATION* material = (COROTATION*)(_tetMesh->materialCopies()[id][materialIndex]);

      MATRIX3 firstPK = material->firstPiolaKirchhoff(_Rs[tetID], _Ss[tetID]);
      forces[0] = firstPK * b[0];
      forces[1] = firstPK * b[1];
      forces[2] = firstPK * b[2];
      forces[3] = firstPK * b[3];

      // add forces to corresponding vertices
      for (int y = 0; y < 4; y++)
      {
        int index = _tetMesh->vertexID(tets[tetID].vertices[y]);
        if (_tetMesh->isSurfaceVertex(index))
        {
          _RCopies[id].segment<3>(index * 3) += forces[y];
        }
      }
    }
  }
  // merge all the copies
  for (int x = 0; x < _tetMesh->totalCores(); x++)
    R += _RCopies[x];

  TIMING_BREAKDOWN::toc("Compute Internal Force");
  return R;
}
