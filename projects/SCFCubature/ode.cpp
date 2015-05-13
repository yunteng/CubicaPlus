#include <iostream>
#include <sstream>
#include <util/IO.h>
#include <util/TIMING_BREAKDOWN.h>
#include <geometry/ODE_BONE.h>
#include <geometry/SKELETON.h>
#include <geometry/SUBSPACE_TET_MESH.h>
#include <geometry/RIGGER.h>
#include <cubature/PENALTY_SCF_CUBATURE_TRAINER.h>
#include <cubature/NNHTP_CUBATURE_GENERATOR.h>


using namespace std;

typedef RIGGER<ODE_BONE> Rigger;
typedef SKELETON<ODE_BONE> Skeleton;

typedef PENALTY_SCF_CUBATURE_TRAINER<ODE_BONE> CubatureApp;
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

    // partition the mesh
    vector<int> tetPartitions;
    rigger->buildSkinningPartition(tetPartitions);
    tetMesh->buildPartitions(tetPartitions);
    tetMesh->loadPartitionBases();

    cubatureApp = new CubatureApp(tetMesh, rigger);
    cubatureGenerator = new CubatureGenerator(cubatureApp);

    // if(partition >= 0 && partition < tetMesh->totalPartitions()){
    for(int partition = 0; partition < tetMesh->totalPartitions(); partition++){

      vector<string> snapshotIdx = IO::getAllSnapshots(dataPath);

      // for each pair of skeletal partitions
      for(int otherPartition = partition + 1; otherPartition < tetMesh->totalPartitions(); otherPartition++){
        cout << "training cubature for partitions " << partition << " and " << otherPartition << endl;
        // set the current training pair
        cubatureApp->setPartitions(partition, otherPartition);

        // group poses with the same relative 
        // skeleton configuration together
        vector<string> trainingPoses;
        cubatureApp->groupPoses(snapshotIdx, trainingPoses);

        // transform the surface position of 
        // otherParition to the local bone 
        // coordinate system of partition
        vector<VECTOR> transformedPositions;
        for(unsigned int x = 0; x < trainingPoses.size(); x++){
          transformedPositions.push_back(cubatureApp->getTransformedColumn(trainingPoses[x]));
        }
        // compute a low rank basis out of 
        // transformedPositions
        MATRIX eigenVectors;
        VECTOR eigenValues;
        MATRIX_UTIL::pcaEigen(transformedPositions, false, eigenVectors, eigenValues);
        int rank = eigenVectors.cols() < 6 ? eigenVectors.cols() : 6;
        eigenVectors.conservativeResize(eigenVectors.rows(), rank);
        MATRIX_UTIL::orthogonalize(eigenVectors);
        MATRIX basisCoordinates(rank, transformedPositions.size());
        for(unsigned int x = 0; x < transformedPositions.size(); x++)
          basisCoordinates.col(x) = eigenVectors.transpose() * transformedPositions[x];

        stringstream cubatureFile("");  
        cubatureFile << tetMesh->filename() << ".scfcubature." << partition << "." << otherPartition;

        // write the basis
        FILE* file = fopen(cubatureFile.str().c_str(), "wb");
        IO::write(eigenVectors, file);
        IO::write(basisCoordinates, file);

        // for each training frame, training the
        // self-collision cubatures
        for(unsigned int x = 0; x < trainingPoses.size(); x++){
          cubatureApp->clearSamples();
          cout << "trainingPose " << trainingPoses[x] << endl;
          
          if(cubatureApp->gatherTrainingData(trainingPoses[x]) == 0){
            int cubatureSize = 0;
            fwrite((void*)&cubatureSize, sizeof(int), 1, file);
            continue;
          }
          // ignore the output from cubatureGenerator
          cubatureApp->setCubatureFilename("./ignore");
          cubatureGenerator->generateCubatures();
          vector<int>& keyPointIDs = cubatureGenerator->keyPointIDs();
          vector<Real>& keyWeights = cubatureGenerator->keyWeights();
          // save the correct cubature format 
          // these cubatures will be loaded from
          // SCF_CUBATURE_LOADER 
          cubatureApp->writePairwiseCubatures(keyPointIDs, keyWeights, file);
        }
        fclose(file);
      }
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
