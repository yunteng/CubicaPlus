/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef SUBSPACE_TET_MESH_H
#define SUBSPACE_TET_MESH_H

#include <SETTINGS.h>
#include <geometry/TET_MESH.h>

using namespace::std;

class SUBSPACE_TET_MESH : public TET_MESH
{
public:
  SUBSPACE_TET_MESH();
  ~SUBSPACE_TET_MESH();

  Real rank() const { return _UBasis.cols(); };
  bool isPartitionedBasis() const { return _usePartitionedBasis; };
  bool isSkinningBasis() const { return _isSkinningBasis; };
  string basisFilename()  const { return _basisFilename; };

  void computePCABasis(Real varianceCutoff);

  MATRIX& U() { return _UBasis; };
  VECTOR& q() { return _q; };

  vector<int>& keyTets()         { return _keyTets; };
  vector<Real>& keyWeights()     { return _keyTetWeights; };

  void updateSubspaceMesh();

  bool readCubatures();

  void generateKeyTetsF();

  VECTOR& reducedInternalForce() { return _reducedR; };
  MATRIX& reducedStiffnessMatrix()    { return _reducedStiffness; };

  MATRIX vertexBasis(int vertexID) 
  { 
    if(vertexID >= _unconstrainedSize){
      MATRIX ret(3, _UBasis.cols());
      ret.setZero();
      return ret;
    }
    return _UBasis.block(vertexID * 3, 0, 3, _UBasis.cols());
  }
  MATRIX vertexBasis(VEC3F* vertex) 
  { 
    int vertexID = _vertexID[vertex];
    return vertexBasis(vertexID);
  }
  MATRIX tetBasis(int tetID)
  {
    MATRIX basis(12, _UBasis.cols());
    for(int x = 0; x < 4; x++){
      basis.block(x * 3, 0, 3, _UBasis.cols()) = vertexBasis(_tets[tetID].vertices[x]);
    }
    return basis;
  }

  void drawKeyTets();

  void computePartitionBases(Real varianceCutoff);
  void loadPartitionBases();
  bool loadPartitionCubatures();
  inline MATRIX& partitionBasis(int x) { return _partitionBases[x]; }

  MATRIX partitionVertexBasis(int partition, int partitionedVertexID)
  {
    return _partitionBases[partition].block(partitionedVertexID * 3, 0, 3, _partitionBases[partition].cols());
  }
  inline int totalPartitionRank() { return _totalPartitionRank; }
  inline int partitionRank(int x)              { assert(x >= 0 && x < _partitionBases.size()); return _partitionBases[x].cols(); }
  inline vector<int>& partitionedKeyTets(int partition) { return _partitionedKeyTets[partition]; };
  inline vector<Real>& partitionedKeyWeights(int partition) { return _partitionedKeyWeights[partition]; };
  inline int partitionRankStartIdx(int x)      { return _partitionRankStartIdx[x]; }
  inline int partitionKeyTetStartIdx(int x) { return _partitionKeyTetStartIdx[x]; }; 
  inline map<pair<int, int>, vector<MATRIX> >& interfaceUs() { return _interfaceUs;}
  inline map<pair<int, int>, vector<MATRIX> >& internalInterfaceUs() { return _internalInterfaceUs;}
  void computeInterfaceUs();
  void computeInternalInterfaceUs();

protected:
  string _basisFilename;
  MATRIX _UBasis;
  bool _isSkinningBasis;
  bool _usePartitionedBasis;
  VECTOR _q;
  VECTOR _reducedR;
  MATRIX _reducedStiffness;

  vector<int> _keyTets;
  vector<Real> _keyTetWeights;

  vector<MATRIX> _partitionBases;
  int _totalPartitionRank;
  vector<int> _partitionRankStartIdx;
  vector<int> _partitionKeyTetStartIdx;
  vector<vector<int> > _partitionedKeyTets;
  vector<vector<Real> > _partitionedKeyWeights;
  map<pair<int, int>, vector<MATRIX> > _interfaceUs;
  map<pair<int, int>, vector<MATRIX> > _internalInterfaceUs;
  int _totalPartitionCubatures;
};

#endif
