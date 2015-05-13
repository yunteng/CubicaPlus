#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
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
    rigger(NULL),
    drawSkeleton(false),
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

    snapshotIds = IO::getAllSnapshots(dataPath);

    startFrame = 0;
    endFrame = snapshotIds.size() - 1;

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

    rigger->readBoneWeights(outputPath + tetmeshName + ".diffusionSkinningWeights");

    rigger->inverseTransformTrainingSamples();

    string transformedDispFilename = dataPath + "restspacedisp.matrix";
    IO::read(transformedDisp, transformedDispFilename);
  }

  void display(){
    rigger->drawBoneSkinning();
    // tetMesh->drawSurfaceFaces();
    if(drawSkeleton)
      skeleton->drawBones();
  }

  void step(){
    if(currentFrame >= 0 && currentFrame <= transformedDisp.cols())
    {
      tetMesh->reset();
      tetMesh->x() = transformedDisp.col(currentFrame);
      tetMesh->updateFullMesh();
    }
  }

  void keyboardFunc(unsigned char key)
  {
    switch(key){
      case 'k':{
        drawSkeleton = !drawSkeleton;
        break;
      }
      case 'i':{
        currentFrame++;
        if(currentFrame > endFrame)
          currentFrame = endFrame;

        cout << "current frame " << snapshotIds[currentFrame] << endl;
        break;
      }
      case 'j':{
        currentFrame--;
        if(currentFrame < 0)
          currentFrame = 0;

        cout << "current frame " << snapshotIds[currentFrame] << endl;
        break;
      }
      case 't':{
        cout << "skinning solution of frame " << snapshotIds[currentFrame] << endl;
        string filename = posePath + skeletonPrefix + snapshotIds[currentFrame] + ".skeleton";
        skeleton->loadFrame(filename);

        bool fromRest = true;

        rigger->updateSkinning(fromRest);
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

  MATRIX transformedDisp;

  vector<string> snapshotIds;
  
  int startFrame;
  int endFrame;
  int skipFrame;
  int currentFrame;
  int previousFrame;

  bool drawSkeleton;

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
