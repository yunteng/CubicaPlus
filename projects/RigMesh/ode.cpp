#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/TET_MESH.h>
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
    rigger(NULL),
    drawSkeleton(false),
    drawConstrained(false),
    currentFrame(0),
    previousFrame(0)
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
    tetMesh = new TET_MESH();

    skeleton = new Skeleton(posePath + skeletonPrefix + "0000.skeleton");

    rigger = new Rigger(skeleton, tetMesh);
    if(!tetMesh->constrained()){
      // constrain the tets penetrated by the skeleton
      rigger->constrainBoneTets();
      delete rigger;
      delete tetMesh;

      // reload the tetmesh since the vertex orderings are changed
      tetMesh = new TET_MESH();
      rigger = new Rigger(skeleton, tetMesh);
    }
    
    if(!rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights")){
      rigger->buildDiffusionSkinning();
      rigger->writeBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");
    }
  }

  void display(){
    rigger->drawBoneSkinning();
    // tetMesh->drawSurfaceFaces();
    if(drawConstrained)
      tetMesh->drawConstrainedNodes();

    if(drawSkeleton)
      skeleton->drawBones();
  }

  void step(){
    if(currentFrame > endFrame)
      return;

    TIMING_BREAKDOWN::startFrame();

    string filename = posePath + skeletonPrefix + IO::itoPaddedString(currentFrame) + ".skeleton";
    if(skeleton->loadFrame(filename)){
      bool fromRest = true;
      TIMING_BREAKDOWN::tic();
      rigger->updateSkinning(fromRest);
      TIMING_BREAKDOWN::toc("Update Skinning");
    }
    
    TIMING_BREAKDOWN::endFrame();

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
        drawConstrained = !drawConstrained;
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
  TET_MESH* tetMesh;
  Skeleton* skeleton;
  Rigger* rigger;

  int startFrame;
  int endFrame;
  int skipFrame;
  int currentFrame;
  int previousFrame;

  bool drawSkeleton;
  bool drawConstrained;

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
