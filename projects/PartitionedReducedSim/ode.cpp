#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/RIGGER.h>
#include <material/PARTITIONED_SUBSPACE_COROTATION_CACHE.h>
#include <integrator/PARTITIONED_SUBSPACE_INTEGRATOR.h>
#include <linearalgebra/NEWTON_RAPHSON.h>
#include <linearalgebra/LU_SOLVER.h>

using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;
typedef PARTITIONED_SUBSPACE_INTEGRATOR<PARTITIONED_SUBSPACE_COROTATION_CACHE, ODE_BONE> Integrator;
typedef NEWTON_RAPHSON<Integrator, LU_SOLVER> Optimizer;

class APPLICATION
{
public:
  APPLICATION(string config):
    configName(config),
    tetMesh(NULL),
    skeleton(NULL),
    rigger(NULL),
    integrator(NULL),
    optimizer(NULL),
    drawSkeleton(false),
    drawSCD(false),
    currentFrame(0),
    previousFrame(0),
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

    simulate   = SIMPLE_PARSER::getBool("simulate subspace", false);

  };
  ~APPLICATION()
  {
    if(optimizer)
      delete optimizer;

    if(integrator)
      delete integrator;

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
    tetMesh->loadPartitionCubatures();

    if(simulate){
      integrator = new Integrator(tetMesh, rigger);
      optimizer = new Optimizer(integrator, &(integrator->hessianInv()));
    }
  }

  void display(){
    rigger->drawBoneSkinning();
    // tetMesh->drawSurfaceFaces();
    
    if(drawSkeleton)
      skeleton->drawBones();
    if(integrator != NULL && drawSCD)
      integrator->scd()->drawSelfCollisionPoints();
  }

  void step(){
    static int stopFrame = SIMPLE_PARSER::getInt("stop frame", endFrame);

    if(currentFrame > endFrame)
      return;

    TIMING_BREAKDOWN::startFrame();

    string filename = posePath + skeletonPrefix + IO::itoPaddedString(currentFrame) + ".skeleton";
    if(skeleton->loadFrame(filename)){
      bool fromRest = false;

      TIMING_BREAKDOWN::tic();
      rigger->updateSkinning(fromRest);
      TIMING_BREAKDOWN::toc("Update Skinning");
    }
    

    if(simulate){
      optimizer->run();
    }
    TIMING_BREAKDOWN::endFrame();

    if(SIMPLE_PARSER::getBool("save mesh", false)){
      tetMesh->writeObj(dataPath + IO::itoPaddedString(currentFrame) + ".obj");
      tetMesh->writeDisplacementFromRest(dataPath + IO::itoPaddedString(currentFrame) + ".state");
    }

    previousFrame = currentFrame;
    currentFrame += skipFrame;
  }

  void keyboardFunc(unsigned char key)
  {
    switch(key){
      case 'k':{
        drawSkeleton = !drawSkeleton;
        break;
      }
      case 'e':{
        drawSCD = !drawSCD;
        break;
      }
      case 'i':{
        currentFrame++;
        break;
      }
      case 'j':{
        currentFrame--;
        break;
      }
      case 'p':{
        TIMING_BREAKDOWN::printTimingBreakdown();
        break;
      }
    }
  }

  void click(VEC3F point)
  {
    
  }

public:
  SUBSPACE_TET_MESH* tetMesh;
  Skeleton* skeleton;
  Rigger* rigger;
  Integrator* integrator;
  Optimizer* optimizer;
  
  int startFrame;
  int endFrame;
  int skipFrame;
  int currentFrame;
  int previousFrame;

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
  string configName(argv[1]);

  VIEWER<APPLICATION>::simulator = new APPLICATION(configName);

  VIEWER<APPLICATION>::simulator->init();
  
  VIEWER<APPLICATION>::init();

  return 0;
}
