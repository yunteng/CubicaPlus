#include <iostream>
#include <sstream>
#include <util/VIEWER.h>
#include <geometry/TET_MESH.h>

using namespace std;

class APPLICATION
{
public:
  APPLICATION(string config):
    configName(config),
    currentFrame(0),
    lastFrame(0)
  {
    if(!SIMPLE_PARSER::parse(configName))
        exit(0);

    outputPath       = SIMPLE_PARSER::getString("output path", "");
    renderPath       = SIMPLE_PARSER::getString("render path", "");
    dataPath         = SIMPLE_PARSER::getString("data path", "");
    posePath         = SIMPLE_PARSER::getString("pose path", "");

    // create working directories
    string mkdirRender = string("mkdir ") + renderPath;
    system(mkdirRender.c_str());
    string mkdirData = string("mkdir ") + dataPath;
    system(mkdirData.c_str());
  };
  ~APPLICATION()
  {

  }
  void init()
  {
    tetMesh = new TET_MESH();
    vector<TET>& tets = tetMesh->tets();
    Real maxVol = -1;
    Real minVol = 1e9;
    Real sumVol = 0;
    for(unsigned int x = 0; x < tets.size(); x++){
      Real vol = tets[x].volume();
      if(vol > maxVol)
        maxVol = vol;
      if(vol < minVol)
        minVol = vol;
      sumVol += vol;
    }
    cout << "======================" << endl;
    cout << "max tet volume " << maxVol << endl;
    cout << "min tet volume " << minVol << endl;
    cout << "avg volume " << sumVol / tets.size() << endl;
    cout << "======================" << endl;
  }

  void display(){
    tetMesh->drawSurfaceFaces();
  }

  void step(){

  }

  void keyboardFunc(unsigned char key)
  {

  }

  void click(VEC3F point)
  {
    
  }

public:
  TET_MESH* tetMesh;

  int currentFrame;
  int lastFrame;

  string configName;

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
