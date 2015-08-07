#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/RIGGER.h>
#include <material/FULLSPACE_COROTATION_CACHE.h>
#include <cubature/PARTITIONED_TET_MESH_CUBATURE_TRAINER.h>
#include <cubature/NNHTP_CUBATURE_GENERATOR.h>


using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;

typedef PARTITIONED_TET_MESH_CUBATURE_TRAINER<FULLSPACE_COROTATION_CACHE, ODE_BONE> CubatureApp;
typedef NNHTP_CUBATURE_GENERATOR<CubatureApp> CubatureGenerator;

class APPLICATION
{
public:
  APPLICATION(string config):
    configName(config),
    tetMesh(NULL),
    skeleton(NULL),
    rigger(NULL),
    cubatureApp(NULL),
    cubatureGenerator(NULL)
  {
    if(!SIMPLE_PARSER::parse(configName))
        exit(0);

    tetmeshName      = SIMPLE_PARSER::getString("tet mesh name", "");
    skeletonPrefix   = SIMPLE_PARSER::getString("skeleton prefix", "");

    outputPath       = SIMPLE_PARSER::getString("output path", "");
    posePath         = SIMPLE_PARSER::getString("pose path", "");

  };
  ~APPLICATION()
  {
    if(cubatureGenerator)
      delete cubatureGenerator;

    if(cubatureApp)
      delete cubatureApp;

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

    rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");

    vector<int> tetPartitions;
    rigger->buildSkinningPartition(tetPartitions);
    tetMesh->buildPartitions(tetPartitions);
    tetMesh->loadPartitionBases();

    cubatureApp = new CubatureApp(tetMesh, rigger);
    cubatureGenerator = new CubatureGenerator(cubatureApp);

    // if(partition >= 0 && partition < tetMesh->totalPartitions()){
    for(int partition = 0; partition < tetMesh->totalPartitions(); partition++){
      cout << "training cubature for partition " << partition << endl;
      cubatureApp->setCurrentPartition(partition);
      cubatureApp->gatherTrainingData();
      string cubatureFilename = tetMesh->filename();
    
      cubatureFilename += ".transformed.partitionedbasis.partition." + IO::intToString(partition) + ".cubature";
      cubatureApp->setCubatureFilename(cubatureFilename);
      cubatureGenerator->generateCubatures();
    }
  }

public:
  SUBSPACE_TET_MESH* tetMesh;
  Skeleton* skeleton;
  Rigger* rigger;

  CubatureApp* cubatureApp;
  CubatureGenerator* cubatureGenerator;

  string configName;

  string tetmeshName;
  string skeletonPrefix;
  string outputPath;
  string posePath;
};

int main(int argc, char* argv[])
{
  if(argc < 2){
    cout << "./bin/PartitionedInternalForceCubature *.cfg" << endl;
    exit(0);
  }

  string configName(argv[1]);

  APPLICATION simulator(configName);

  simulator.init();

  return 0;
}
