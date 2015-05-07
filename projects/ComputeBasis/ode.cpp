#include <iostream>
#include <sstream>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/RIGGER.h>

using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;

class APPLICATION
{
public:
  APPLICATION(string config):
    configName(config),
    tetMesh(NULL),
    skeleton(NULL),
    rigger(NULL)
  {
    if(!SIMPLE_PARSER::parse(configName))
        exit(0);

    tetmeshName      = SIMPLE_PARSER::getString("tet mesh name", "");
    skeletonPrefix   = SIMPLE_PARSER::getString("skeleton prefix", "");

    outputPath       = SIMPLE_PARSER::getString("output path", "");
    dataPath         = SIMPLE_PARSER::getString("data path", "");
    posePath         = SIMPLE_PARSER::getString("pose path", "");

  };
  ~APPLICATION()
  {
    if(rigger)
      delete rigger;

    if(skeleton)
      delete skeleton;

    if(tetMesh)
      delete tetMesh;

  }
  void init()
  {
    tetMesh = new SUBSPACE_TET_MESH();

    skeleton = new Skeleton(posePath + skeletonPrefix + "0000.skeleton");

    rigger = new Rigger(skeleton, tetMesh);

    // rigger->buildRigidSkinning();
    rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");

    Real variance = SIMPLE_PARSER::getFloat("pca variance", 0.9999);
    bool usePartitionedBasis = SIMPLE_PARSER::getBool("partitioned basis", false);
    
    if(usePartitionedBasis){
      vector<int> tetPartitions;
      rigger->buildSkinningPartition(tetPartitions);
      tetMesh->buildPartitions(tetPartitions);
      tetMesh->computePartitionBases(variance);  
    }else{
      tetMesh->computePCABasis(variance);
    }
  }

public:
  string configName;

  SUBSPACE_TET_MESH* tetMesh;
  Skeleton* skeleton;
  Rigger* rigger;

  string tetmeshName;
  string skeletonPrefix;
  string outputPath;
  string dataPath;
  string posePath;
};

int main(int argc, char* argv[])
{
  string configName(argv[1]);

  APPLICATION simulator(configName);

  simulator.init();

  return 0;
}
