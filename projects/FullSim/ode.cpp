#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/TET_MESH.h>
#include <geometry/RIGGER.h>
#include <material/FULLSPACE_COROTATION_CACHE.h>
#include <integrator/FULLSPACE_INTEGRATOR.h>
#include <linearalgebra/NEWTON_RAPHSON.h>
#include <linearalgebra/JACOBI_PRECONDITIONER.h>
#include <linearalgebra/CONJUGATE_GRADIENT.h>


using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;
typedef FULLSPACE_INTEGRATOR<FULLSPACE_COROTATION_CACHE, ODE_BONE> Integrator;
typedef CONJUGATE_GRADIENT<JACOBI_PRECONDITIONER, Integrator> CGSolver;
typedef NEWTON_RAPHSON<Integrator, CGSolver> Optimizer;

class APPLICATION
{
public:
  APPLICATION(string config):
    configName(config),
    tetMesh(NULL),
    skeleton(NULL),
    rigger(NULL),
    integrator(NULL),
    preconditioner(NULL),
    cgSolver(NULL),
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

    simulate   = SIMPLE_PARSER::getBool("simulate full", false);

  };
  ~APPLICATION()
  {
    if(optimizer)
      delete optimizer;

    if(preconditioner)
      delete preconditioner;

    if(cgSolver)
      delete cgSolver;

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
    tetMesh = new TET_MESH();

    skeleton = new Skeleton(posePath + skeletonPrefix + "0000.skeleton");

    rigger = new Rigger(skeleton, tetMesh);

    rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");

    if(simulate){
      integrator = new Integrator(tetMesh, rigger);
      preconditioner = new JACOBI_PRECONDITIONER(integrator->systemMatrixDiag());
      cgSolver = new CGSolver(integrator, preconditioner);
      optimizer = new Optimizer(integrator, cgSolver);
    }

    // let the simulation start from a non-rest pose
    int preloadFrame = SIMPLE_PARSER::getInt("preload frame", 0);
    if(preloadFrame != 0){
      tetMesh->reset();
      string filename = posePath + skeletonPrefix + IO::itoPaddedString(preloadFrame) + ".skeleton";
      skeleton->loadFrame(filename);
      skeleton->fixSkeletonStructure();
      Real fromRest = true;
      // just in case there is no precompute displacement for this frame
      // update using skinning weights first
      rigger->updateSkinning(fromRest);

      tetMesh->readDisplacementFromRest(dataPath + IO::itoPaddedString(preloadFrame) + ".state");
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
    if(currentFrame > endFrame)
      exit(0);

    TIMING_BREAKDOWN::startFrame();

    // load skeleton
    string filename = posePath + skeletonPrefix + IO::itoPaddedString(currentFrame) + ".skeleton";
    if(skeleton->loadFrame(filename)){
      bool fromRest = false;

      TIMING_BREAKDOWN::tic();
      rigger->updateSkinning(fromRest);
      TIMING_BREAKDOWN::toc("Update Skinning");
    }
    

    if(simulate){
      optimizer->run();
    }else{
      tetMesh->readDisplacementFromRest(dataPath + IO::itoPaddedString(currentFrame) + ".state");
    }
    
    TIMING_BREAKDOWN::endFrame();

    if(simulate){
      tetMesh->writeDisplacementFromRest(dataPath + IO::itoPaddedString(currentFrame) + ".state");
      integrator->writeSelfCollisionResponses(dataPath + IO::itoPaddedString(currentFrame) + ".collisionresponse");
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
  Integrator* integrator;
  JACOBI_PRECONDITIONER* preconditioner;
  CGSolver* cgSolver;
  Optimizer* optimizer;
  
  int startFrame;
  int endFrame;
  int previousFrame;
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
bool VIEWER<T>::animate = true;

template <class T>
bool VIEWER<T>::step = false;

int main(int argc, char* argv[])
{
  string configName(argv[1]);

  VIEWER<APPLICATION>::simulator = new APPLICATION(configName);

  VIEWER<APPLICATION>::simulator->init();
  
  VIEWER<APPLICATION>::init();

  return 0;
}
