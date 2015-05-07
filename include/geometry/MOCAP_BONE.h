#ifndef MOCAP_BONE_H
#define MOCAP_BONE_H

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <SETTINGS.h>
#include <vector>
#include <util/SIMPLE_PARSER.h>
#include <geometry/DUAL_QUATERNION.h>

using namespace::std;

class MOCAP_BONE
{
public:
  // MOCAP_BONE(MOCAP_BONE* parent, const MATRIX4& transform, Real length = 1.0){

  // };
  // MOCAP_BONE(int parentId, const MATRIX4& transform, Real length = 1.0){

  // };
  // MOCAP_BONE():
  //   _parent(NULL),
  //   _parentId(-1) 
  // {
  // };
  
  MOCAP_BONE(FILE* file):
    _parent(NULL),
    _parentId(-1)
  {
    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);

    float trans[3];

    fscanf(file,"%i %f %f %f\n", &_id, &trans[0], &trans[1], &trans[2]);
    _restEndPoint = VEC3F(trans[0], trans[1], trans[2]) * scale;

    fscanf(file,"%f %f %f\n", &trans[0], &trans[1], &trans[2]);
    _restBeginPoint = VEC3F(trans[0], trans[1], trans[2]) * scale;

    _restTranslation.setZero();
    _previousTranslation.setZero();

    _restRotation.setIdentity();
    _restRotationInv.setIdentity();
    _previousRotation.setIdentity();

    _translation.setZero();
    _rotation.setIdentity();

    _length = (_restEndPoint - _restBeginPoint).norm();

    _dqFromRest = DUAL_QUATERNION(_rotation, _translation);
  };
  void write(FILE* file)
  {
    /*Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);

    float trans[3];
    trans[0] = _translation[0] / scale;
    trans[1] = _translation[1] / scale;
    trans[2] = _translation[2] / scale;

    fprintf(file,"%i %f %f %f %i\n", _id, trans[0], trans[1], trans[2], _parentId);

    float quat[4];
    quat[0] = _rotation.x();
    quat[1] = _rotation.y();
    quat[2] = _rotation.z();
    quat[3] = _rotation.w();

    fprintf(file,"%f %f %f %f\n", quat[0], quat[1], quat[2], quat[3]);

    float length = _length / scale;
    float radius = _radius / scale;
    fprintf(file,"%f %f\n", length, radius);*/
  }
  void updateTransform(FILE* file)
  {
    _previousTranslation = _translation;
    _previousRotation = _rotation;

    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);
    int index;
    float trans[3];
    fscanf(file,"%i %f %f %f\n", &index, &trans[0], &trans[1], &trans[2]);
    _translation = VEC3F(trans[0], trans[1], trans[2]) * scale;

    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);

    _rotation = QUATERNION(quat[0], quat[1], quat[2], quat[3]);
    _rotation.normalize();

    _dqFromRest = DUAL_QUATERNION(_rotation, _translation);
  }

  void updateTransformAccordingToParent()
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }

  void updateAccordingToParent()
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }

  void interpolateFromRestToRecord(Real step, int recordID)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }

  void updateRecentRecordTransformAccordingToParent()
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }
  
  inline void setParentID(int id) { _parentId = id; };
  inline void setParent(MOCAP_BONE* bone) { _parent = bone; }; 
  
  inline VEC3F& translation() { return _translation; };
  inline MATRIX3 rotation() { return _rotation.toRotationMatrix(); };

  inline VEC3F& restTranslation() { return _restTranslation; };
  inline MATRIX3 restRotation() { return _restRotation.toRotationMatrix(); };

  inline MATRIX3& restRotationInv() { return _restRotationInv; };

  void reset() {
    _translation = _restTranslation;
    _rotation = _restRotation;
  };

  inline pair<VEC3F, VEC3F> restBoneSegments()
  {
    return make_pair(_restEndPoint, _restBeginPoint);
  }

  void draw() {
    VEC3F beginVertex = _rotation._transformVector(_restBeginPoint) + _translation;

    VEC3F endVertex = _rotation._transformVector(_restEndPoint) + _translation;

    VEC3F centerVertex = (beginVertex + endVertex) * 0.5;

    glPointSize(10.0f);
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
      glVertex3f(beginVertex[0], beginVertex[1], beginVertex[2]);
    glColor3f(0.0, 0.0, 1.0);
      glVertex3f(endVertex[0], endVertex[1], endVertex[2]);
    glEnd();

    glLineWidth(3.0f);
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
      glVertex3f(beginVertex[0], beginVertex[1], beginVertex[2]);
      glVertex3f(endVertex[0], endVertex[1], endVertex[2]);
    glEnd();

    glPointSize(10.0f);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
      glVertex3f(centerVertex[0], centerVertex[1], centerVertex[2]);
    glEnd();
  };

  inline VEC3F transformFromRest(const VEC3F& pos){
    return _rotation._transformVector(pos) + _translation;
  }

  inline VEC3F transformFromPrevious(const VEC3F& pos){
    VEC3F ret = _previousRotation.conjugate()._transformVector(pos - _previousTranslation);
    ret = _rotation._transformVector(ret) + _translation;
    return ret;
  }

  inline void transformFromRestUsingRecords(const VEC3F& pos, vector<VEC3F>& output)
  {
    output.clear();
    for(unsigned int x = 0; x < _recordedTranslations.size(); x++){
      output.push_back(_recordedRotations[x].toRotationMatrix() * pos + _recordedTranslations[x]);
    }
  }
  
  inline DUAL_QUATERNION dualQuaternionRecordFromRest(int recordID)
  {
    assert(recordID >= 0 && recordID < _recordedRotations.size());
    
    return DUAL_QUATERNION(_recordedRotations[recordID], _recordedTranslations[recordID]);
  }

  inline DUAL_QUATERNION dualQuaternionFromRest()
  {
    /*QUATERNION relativeRotation = _rotation * _restRotation.conjugate();

    VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_restTranslation);
    return DUAL_QUATERNION(relativeRotation, relativeTranslation);*/

    // return DUAL_QUATERNION(_rotation, _translation);
    return _dqFromRest;
  }

  inline DUAL_QUATERNION dualQuaternionFromPrevious()
  {
    QUATERNION relativeRotation = _rotation * _previousRotation.conjugate();

    VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_previousTranslation);
    return DUAL_QUATERNION(relativeRotation, relativeTranslation);
  }

  inline MATRIX3 getRotationRecord(int recordID){
    assert(recordID >= 0 && recordID < _recordedRotations.size());
    return _recordedRotations[recordID].toRotationMatrix();
  }
  inline VEC3F& getTranslationRecord(int recordID){
    assert(recordID >= 0 && recordID < _recordedTranslations.size());
    return _recordedTranslations[recordID];
  }
  void recordTransform(FILE* file)
  {
    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);
    int index;
    float trans[3];

    fscanf(file,"%i %f %f %f\n", &index, &trans[0], &trans[1], &trans[2]);

    QUATERNION rotation;
    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
    rotation.w() = quat[0]; rotation.x() = quat[1]; 
    rotation.y() = quat[2]; rotation.z() = quat[3];
    rotation.normalize();

    _recordedTranslations.push_back(VEC3F(trans[0], trans[1], trans[2]) * scale);
    _recordedRotations.push_back(rotation);
  };

  inline void clearRecords(){
    _recordedTranslations.clear();
    _recordedRotations.clear();
  };

  bool isRestPose(){
    QUATERNION diffRot = _rotation - _restRotation;
    return _translation.norm() + diffRot.x() * diffRot.x() + diffRot.y() * diffRot.y() + diffRot.z() * diffRot.z() < 1e-5;
  }

  bool isRecordPose(int recordID){
    if(recordID < 0 && recordID >= _recordedRotations.size())
      return false;
    QUATERNION diffRot = _rotation - _recordedRotations[recordID];
    return (_translation - _recordedTranslations[recordID]).norm() + diffRot.x() * diffRot.x() + diffRot.y() * diffRot.y() + diffRot.z() * diffRot.z() < 1e-5;
  }

  void computeRelativeRT(MOCAP_BONE* otherBone, VEC3F& relativeTranslation, QUATERNION& relativeRotation)
  {
    relativeRotation = _rotation * (otherBone->_rotation).conjugate();
    if(relativeRotation.w() < 0){
      relativeRotation.x() *= -1;
      relativeRotation.y() *= -1;
      relativeRotation.z() *= -1;
      relativeRotation.w() *= -1;
    }
    relativeTranslation = _rotation.conjugate()._transformVector(_translation - otherBone->translation());
  }


private:
  VEC3F _restBeginPoint;
  VEC3F _restEndPoint;

  VEC3F _translation;
  QUATERNION _rotation;

  VEC3F _previousTranslation;
  QUATERNION _previousRotation;

  DUAL_QUATERNION _dqFromRest;

  vector<VEC3F>       _recordedTranslations;
  vector<QUATERNION>  _recordedRotations;

  VEC3F _restTranslation;
  QUATERNION _restRotation;

  MATRIX3 _restRotationInv;

  Real _length;

  MOCAP_BONE* _parent;
  int _id;
  int _parentId;
};

#endif
