#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/SPHERE.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/RIGGER.h>
#include <material/PARTITIONED_FULLSPACE_COROTATION_CACHE.h>
#include <material/PARTITIONED_SUBSPACE_COROTATION_CACHE.h>
#include <integrator/PARTITIONED_HYBRID_INTEGRATOR.h>
#include <linearalgebra/HYBRID_NEWTON_RAPHSON.h>
#include <linearalgebra/HYBRID_PRECONDITIONER.h>
#include <linearalgebra/CONJUGATE_GRADIENT.h>

using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;
typedef PARTITIONED_HYBRID_INTEGRATOR<PARTITIONED_FULLSPACE_COROTATION_CACHE, PARTITIONED_SUBSPACE_COROTATION_CACHE, ODE_BONE> Integrator;
typedef HYBRID_PRECONDITIONER<Integrator> Preconditioner;
typedef CONJUGATE_GRADIENT<Preconditioner, Integrator> CGSolver;
typedef HYBRID_NEWTON_RAPHSON<Integrator, CGSolver> Optimizer;

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
    drawFullsimRegion(false),
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
    string mkdirRender = string("mkdir -p ") + renderPath;
    system(mkdirRender.c_str());
    string mkdirData = string("mkdir -p ") + dataPath;
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
    tetMesh = new SUBSPACE_TET_MESH();
    if(SIMPLE_PARSER::getBool("save mesh", false)){
      tetMesh->writeObj(renderPath + "0000.obj");
    }

    skeleton = new Skeleton(posePath + skeletonPrefix + "0000.skeleton");

    rigger = new Rigger(skeleton, tetMesh);
    
    /*
    read skinning weights
    */
    rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");

    /*
    Partition the mesh, associate each tet with the bone that has maximum weights on its vertices
    */
    vector<int> tetPartitions;
    rigger->buildSkinningPartition(tetPartitions);
    tetMesh->buildPartitions(tetPartitions);
    /*
    Load basis for each partition
    */
    tetMesh->loadPartitionBases();
    /*
    Load material force cubatures for each 
    partition, the self collision cubatures
    are loaded in self collision detector
    */
    tetMesh->loadPartitionCubatures();

    if(simulate){
      integrator = new Integrator(tetMesh, rigger);
      /*
      initialize PCG
      */
      preconditioner = new Preconditioner(integrator);
      cgSolver = new CGSolver(integrator, preconditioner);
      /*
      the Newton solver
      */
      optimizer = new Optimizer(integrator, cgSolver);
    }
    /*
    allow the simulation to start from a precomputed non-rest pose
    */
    int preloadFrame = SIMPLE_PARSER::getInt("preload frame", 0);
    if(preloadFrame != 0){
      tetMesh->reset();
      string filename = posePath + skeletonPrefix + IO::itoPaddedString(preloadFrame) + ".skeleton";
      skeleton->loadFrame(filename);
      Real fromRest = true;
      // just in case there is no precompute displacement for this frame
      // update using skinning weights first
      rigger->updateSkinning(fromRest);

      tetMesh->readDisplacementFromRest(dataPath + IO::itoPaddedString(preloadFrame) + ".state");
    }
    cout << "total dofs " << tetMesh->partitionDofStartIdx(tetMesh->totalPartitions()) << endl;
  }

  void display(){
    rigger->drawBoneSkinning();
    // tetMesh->drawSurfaceFaces();

    if(drawFullsimRegion)
      // skeleton->drawBones();
      tetMesh->drawFullsimVertices();

    // draw self collision vertices
    if(integrator != NULL && drawSCD){
      integrator->scd()->drawSelfCollisionPoints();
    }
  }

  void step(){
    if(currentFrame > endFrame)
      return;

    TIMING_BREAKDOWN::startFrame();

    // load skeleton and perform skinning update
    string filename = posePath + skeletonPrefix + IO::itoPaddedString(currentFrame) + ".skeleton";
    if(skeleton->loadFrame(filename)){

      TIMING_BREAKDOWN::tic();
      bool fromRest = false;
      rigger->updateSkinning(fromRest);
      TIMING_BREAKDOWN::toc("Update Skinning");
    }

    if(simulate){
      optimizer->run();
    }
    
    TIMING_BREAKDOWN::endFrame();

    if(SIMPLE_PARSER::getBool("save mesh", false)){
      tetMesh->writeObj(renderPath + IO::itoPaddedString(currentFrame) + ".obj");
      tetMesh->writeDisplacementFromRest(dataPath + IO::itoPaddedString(currentFrame) + ".state");
    }
      
    previousFrame = currentFrame;
    currentFrame += skipFrame;
  }

  void keyboardFunc(unsigned char key)
  {
    switch(key){
      case 'k':{
        drawFullsimRegion = !drawFullsimRegion;
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
  Preconditioner* preconditioner;
  CGSolver* cgSolver;
  Optimizer* optimizer;
  
  int startFrame;
  int endFrame;
  int skipFrame;
  int currentFrame;
  int previousFrame;

  bool simulate;

  bool drawFullsimRegion;
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
