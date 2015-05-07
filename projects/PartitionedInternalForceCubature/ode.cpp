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
    cubatureGenerator(NULL),
    drawSkeleton(false),
    drawSCD(false),
    currentFrame(0),
    simulate(false)
  {
    if(!SIMPLE_PARSER::parse(configName))
        exit(0);

    tetmeshName      = SIMPLE_PARSER::getString("tet mesh name", "");
    skeletonPrefix   = SIMPLE_PARSER::getString("skeleton prefix", "");

    outputPath       = SIMPLE_PARSER::getString("output path", "");
    renderPath       = SIMPLE_PARSER::getString("render path", "");
    dataPath         = SIMPLE_PARSER::getString("data path", "");
    posePath         = SIMPLE_PARSER::getString("pose path", "");

    // create working directories
    string mkdirRender = string("mkdir ") + renderPath;
    system(mkdirRender.c_str());
    string mkdirData = string("mkdir ") + dataPath;
    system(mkdirData.c_str());

    startFrame = SIMPLE_PARSER::getInt("start frame", 1);
    endFrame   = SIMPLE_PARSER::getInt("snapshots", 100);
    skipFrame  = SIMPLE_PARSER::getInt("skip frame", 1);
    currentFrame = startFrame;

    simulate   = SIMPLE_PARSER::getBool("simulate full", false);

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
  
  int startFrame;
  int endFrame;
  int skipFrame;
  int currentFrame;

  bool simulate;

  bool drawSkeleton;
  bool drawSCD;

  string configName;

  string tetmeshName;
  string skeletonPrefix;
  string outputPath;
  string renderPath;
  string dataPath;
  string posePath;
};


// declare static member variables
template <class T> 
GLVU VIEWER<T>::glvu;

template <class T> 
T* VIEWER<T>::simulator = NULL;

template <class T>
bool VIEWER<T>::animate = false;

template <class T>
bool VIEWER<T>::step = true;

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
