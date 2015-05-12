#ifndef RIGGER_H
#define RIGGER_H

#include <SETTINGS.h>
#include <geometry/SKELETON.h>
#include <geometry/TET_MESH.h>
#include <util/SIMPLE_PARSER.h>

template<class BONE>
class RIGGER
{
public:
  RIGGER(SKELETON<BONE>* skeleton, TET_MESH* tetMesh);

  /*
  access
  */
  SKELETON<BONE>* skeleton() { return _skeleton; };
  inline vector<vector<pair<int, Real> > >& skinning() { return _skinning; };
  inline string skinningMethod() const { return _skinningMethod; };

  inline VECTOR& skinningDisp() {return _skinningDisp; };
  inline vector<MATRIX3>& skinningRotation() {return _skinningRotation; };
  inline vector<int>& nearestBones() { return _nearestBones; }

  /*
  draw function
  */
  void drawBoneSkinning();
  void drawBoneSkinning(int boneID);
  /*
  rigging functions
  */
  void constrainBoneTets();
  void buildRigidSkinning();
  void buildSkinningPartition(vector<int>& tetPartitions);

  // update the mesh from previous pose or the rest pose
  void updateSkinning(bool fromRest)
  {
    if(_skinningMethod.compare("dual quaternion") == 0)
    {
      updateDualQuaternionSkinning(fromRest);
    }
    else
    {
      cout << "UNKNOWN skinning method " << _skinningMethod << endl;
      exit(0);
    }
  }

  void updateDualQuaternionSkinning(bool fromRest);

  /*
  compute diffusion skinning weights
  */
  void buildDiffusionSkinning();
  
  /*
  partition either the original mesh or its low-res embedding, used for collision detection
  */
  void buildSkinningPartition(vector<vector<int> >& partitionSurfaceVertices, vector<vector<int> >& partitionTets, bool useLowresTets);

  /*
  pull back the training samples to before-skinning space
  */
  void inverseTransformTrainingSamples();

  /*
  IO stuff
  */
  void writeBoneWeights(const string& filename);
  bool readBoneWeights(const string& filename);

private:
  void computeDQblend(int vertexID, DUAL_QUATERNION& dq_blend);
  /*
  diffusion skinning
  */
  Real conductionWeight(Real conductance, Real distToNearest);
  void computeConduction(Real conductance, SpMat& matrix);
  
  void computeConductionRHS(int boneIndex, Real conductance, VECTOR& rhs);

  /*
  normalize the weight sum for each vertex to 1
  */
  void normalizeWeights();

  struct orderByValueGreaterThan{
    bool operator()(const pair<int, Real>& a, const pair<int, Real>& b) { return a.second > b.second; }
  }sortWeight;

private:
  SKELETON<BONE>* _skeleton;
  TET_MESH* _tetMesh;

  vector<int> _nearestBones;
  vector<Real> _distsToNearest;

  vector<DUAL_QUATERNION> _previousDq_blend;

  string _skinningMethod;
  vector<vector<pair<int, Real> > > _skinning;
  vector<int> _maxWeightIndex;

  VECTOR _skinningDisp;
  vector<MATRIX3> _skinningRotation;
};

#include "RIGGER.inl"
#include "RIGGER_DIFFUSE.inl"
#endif
