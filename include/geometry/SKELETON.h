#ifndef SKELETON_H
#define SKELETON_H

#include <SETTINGS.h>
#include <iostream>
#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <util/RGB_HSV.h>

using namespace::std;

template<class BONE>
class SKELETON
{
public:
  SKELETON();
  SKELETON(const string& filename);
  ~SKELETON();

  void loadSturcture(const string& filname);

  bool loadFrame(const string& filename);

  void writeFrame(const string& filename);

  void drawBones();
  void fixSkeletonStructure();
  
  int totalBones()              { return _bones.size();};
  vector<BONE*>& bones()        { return _bones; };
  const vector<VEC3F>& colors() { return _colors; }

  /*
  dealing with weights training
  */
  void recordFrame(const string& filename);
  void clearRecords();
  void getTransformFromRestUsingRecords(int boneID, const VEC3F& pos, vector<VEC3F>& output);

  inline MATRIX3 getBoneRotationRecord(int boneID, int recordID){
    return _bones[boneID]->getRotationRecord(recordID);
  }
  inline VEC3F getBoneTranslationRecord(int boneID, int recordID){
    return _bones[boneID]->getTranslationRecord(recordID);
  }
  bool isRestPose()
  {
    for(unsigned int x = 0; x < _bones.size(); x++)
      if(!_bones[x]->isRestPose())
        return false;
    return true;
  }
  bool isRecordPose(int recordID)
  {
    for(unsigned int x = 0; x < _bones.size(); x++)
      if(!_bones[x]->isRecordPose(recordID))
        return false;
    return true;
  }
  int numberOfRecords() { return _numberOfRecords; };

  void interpolateFromRestToRecord(int step, int total, int recordID);

  void interpolateFromPreviousToCurrentRecord(int step, int total, int recordID);

  void computeRelativeTransform(int x, int y, VEC3F& relativeTranslation, QUATERNION& relativeRotation);
  VEC3F computeRestJointPosition(int x, int y);

private:
  void computeColors();
  void fixRecentRecordStructure();
  
private:
  vector<BONE*> _bones;
  vector<VEC3F> _colors;
  vector<vector<pair<int, Real> > > _skinning;
  vector<int> _boneHierarchy;
  int _numberOfRecords;
};

#include "SKELETON.inl"

#endif
