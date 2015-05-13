template<class BONE>
Real RIGGER<BONE>::conductionWeight(Real conductance, Real distToNearest)
{
  distToNearest = max(distToNearest, 1e-4);
  return conductance / (distToNearest * distToNearest);
}
template<class BONE>
void RIGGER<BONE>::computeConduction(Real conductance, SpMat& matrix)
{
  for (int i = 0; i < _distsToNearest.size(); ++i)
  {
    matrix.coeffRef(i, i) += conductionWeight(conductance, _distsToNearest[i]);
  }
}

template<class BONE>
void RIGGER<BONE>::computeConductionRHS(int boneIndex, Real conductance, VECTOR& rhs)
{
  for(int i = 0; i < _distsToNearest.size(); ++i) {
    if(_nearestBones[i] == boneIndex)
      rhs[i] = conductionWeight(conductance, _distsToNearest[i]);
  }
}

template<class BONE>
void RIGGER<BONE>::buildDiffusionSkinning()
{
  buildRigidSkinning();

  cout << " Compute diffusion skinning weights..." << endl;

  cout << "  compute vertex adjacency..."; flush(cout);
  vector<vector<int> > adjacentVerts;
  _tetMesh->computeVertexAdjacency(adjacentVerts);
  cout << " done." << endl;

  _skinning.clear();
  _skinning.resize(_tetMesh->vertices().size());

  cout << "  compute laplacian matrix..."; flush(cout);
  SpMat systemMatrix;
  _tetMesh->computeLaplacianMatrix(adjacentVerts, systemMatrix);
  cout << " done." << endl;

  cout << "  compute bone conductance..."; flush(cout);
  Real boneConductance = SIMPLE_PARSER::getFloat("bone conductance", 10.0);
  computeConduction(boneConductance, systemMatrix);
  cout << " done." << endl;

  
  systemMatrix.makeCompressed();

  #pragma omp critical
  {
  cout << "  factorize system matrix..."; flush(cout);
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<SpMat::Index> > sparseluSolver;
  sparseluSolver.analyzePattern(systemMatrix);
  sparseluSolver.factorize(systemMatrix);
  cout << " done." << endl;

  
  cout << "  solve individual bone weights..."; flush(cout);
  VECTOR rhs(systemMatrix.cols());
  VECTOR x(systemMatrix.cols());
  for(int i = 0; i < _skeleton->totalBones(); i++){
    cout << i << " "; flush(cout);
    rhs.setZero();
    computeConductionRHS(i, boneConductance, rhs);
    x = sparseluSolver.solve(rhs);
    for(unsigned int j = 0; j < x.size(); j++){
      Real w = x[j];
      if(!(w <= 1.0 + 1e-7 && w >= 0.0 - 1e-7)) {
        cout << "Got invalid weight value (" << j << "): " << w << endl;
      }
      w = std::min(std::max(w, 0.0), 1.0);
      _skinning[j].push_back(make_pair(i, w));
    } 
  }
  cout << " done." << endl;
  }

  cout << "  normalize bone weights..."; flush(cout);
  normalizeWeights();
  cout << " done." << endl;

  writeBoneWeights(SIMPLE_PARSER::getString("output path", "") + SIMPLE_PARSER::getString("tet mesh name", "") + ".diffusionSkinningWeights");

}
