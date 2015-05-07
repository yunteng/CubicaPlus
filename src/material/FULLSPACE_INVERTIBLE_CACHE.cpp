#include <material/FULLSPACE_INVERTIBLE_CACHE.h>
#include <util/MATRIX_UTIL.h>
#include <util/TIMING_BREAKDOWN.h>

FULLSPACE_INVERTIBLE_CACHE::FULLSPACE_INVERTIBLE_CACHE(TET_MESH* tetMesh):
  _tetMesh(tetMesh)
{
  int totalTets = tetMesh->tets().size();
  _Us.resize(totalTets);
  _Vs.resize(totalTets);
  _Fhats.resize(totalTets);
  _stiffnesses.resize(totalTets);

  int rank = tetMesh->dofs();
  _RCopies = new VECTOR[tetMesh->totalCores()];
  for(int x = 0; x < tetMesh->totalCores(); x++)
    _RCopies[x].resize(rank);
}
FULLSPACE_INVERTIBLE_CACHE::~FULLSPACE_INVERTIBLE_CACHE()
{
  delete[] _RCopies;
}

void FULLSPACE_INVERTIBLE_CACHE::cacheDecompositions()
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
      INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);
      
      MATRIX3& F = Fs[x];
      Real* stiffnessData = _stiffnesses[x].data();
      material->diagonalizeF(F, _Us[x], _Fhats[x], _Vs[x]);
      material->stiffnessDensity(_Us[x], _Fhats[x], _Vs[x], stiffnessData);
    }
  }
  TIMING_BREAKDOWN::toc("Cache Material Decompositions");
  
}
void FULLSPACE_INVERTIBLE_CACHE::computeIndividualTetInternalForces(VECTOR& output)
{

}

void FULLSPACE_INVERTIBLE_CACHE::computeIndividualTetInternalForces(vector<int>& tetIDs, VECTOR& output)
{

}

VECTOR& FULLSPACE_INVERTIBLE_CACHE::computeInternalForce()
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
      INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);

      MATRIX3 firstPK = material->firstPiolaKirchhoff(_Us[x], _Fhats[x], _Vs[x]);
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

Real FULLSPACE_INVERTIBLE_CACHE::computeElasticEnergy()
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
      
      INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);

      elasticEnergy += tets[x].restVolume() * material->strainEnergy(tets[x].F());
    }
  }
  TIMING_BREAKDOWN::toc("Compute Elastic Energy");
  return elasticEnergy;
}

COO_MATRIX& FULLSPACE_INVERTIBLE_CACHE::computeStiffnessMatrix()
{
  TIMING_BREAKDOWN::tic();

  COO_MATRIX& stiffness = _tetMesh->stiffnessMatrix();

  if(stiffness.rows() != _tetMesh->dofs())
    stiffness.resize(_tetMesh->dofs(), _tetMesh->dofs());

  vector<TET>& tets = _tetMesh->tets();

  vector<int> offsets(_tetMesh->totalCores() + 1);

#if USING_OPENMP
#pragma omp parallel
#endif
  { 
    #if USING_OPENMP
    const int id  = omp_get_thread_num();
    #else
    const int id = 0;
    #endif
    
    COO_MATRIX localCooStiffness;

    // vector<TRIPLET>& mat = localCooStiffness.matrix();

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
      
      INVERTIBLE* material = (INVERTIBLE*)(_tetMesh->materialCopies()[id][materialIndex]);
      
      // MATRIX3& R = _Rs[x];
      // MATRIX3& S = _Ss[x];
      // MATRIX3& L = _Ls[x];
      MATRIX9& diagonalStiffness = _stiffnesses[x];
      MATRIX3& U = _Us[x];
      MATRIX3& V = _Vs[x];
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

        // rotate deltaF
        MATRIX3 rotated = U.transpose() * MATRIX_UTIL::repackVec9(deltaF) * V;
        deltaF = MATRIX_UTIL::flatternMatrix3(rotated);
        
        VECTOR contraction = diagonalStiffness  * deltaF;
        MATRIX3 deltaP = MATRIX_UTIL::repackVec9(contraction);

        // rotate deltaP back
        deltaP = U * deltaP * V.transpose();

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
            localCooStiffness.add(forceVecs[z][a], row + a, col);
            // mat.push_back(TRIPLET())
              // mat.push_back(make_pair(make_pair(row + a, col), forceVecs[z][a]));
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

