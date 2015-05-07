#ifndef ODE_BONE_H
#define ODE_BONE_H

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

class ODE_BONE
{
public:
  // ODE_BONE(ODE_BONE* parent, const MATRIX4& transform, Real length = 1.0){

  // };
  // ODE_BONE(int parentId, const MATRIX4& transform, Real length = 1.0){

  // };
  // ODE_BONE():
  //   _parent(NULL),
  //   _parentId(-1) 
  // {
  // };
  
  ODE_BONE(FILE* file):
    _parent(NULL),
    _parentId(-1)
  {
    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);
    // int index;
    float trans[3];
    int previousIndex;
    fscanf(file,"%i %f %f %f %i\n", &_id, &trans[0], &trans[1], &trans[2], &previousIndex);

    QUATERNION rotation;
    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
    rotation.x() = quat[0]; rotation.y() = quat[1]; 
    rotation.z() = quat[2]; rotation.w() = quat[3];
    rotation.normalize();
    float length, radius;
    fscanf(file,"%f %f\n", &length, &radius);

    _length = length * scale;
    _radius = radius * scale;

    _restTranslation = _previousTranslation = _translation = VEC3F(trans[0], trans[1], trans[2]) * scale;

    _restRotation = _previousRotation = _rotation = rotation;

    _restRotationInv = _restRotation.toRotationMatrix().transpose();

    _parentId = previousIndex;
    
    QUATERNION relativeRotation = _rotation * _restRotation.conjugate();

    VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_restTranslation);

    _dqFromRest = DUAL_QUATERNION(relativeRotation, relativeTranslation);

  };
  void write(FILE* file)
  {
    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);

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
    fprintf(file,"%f %f\n", length, radius);
  }
  void updateTransform(FILE* file)
  {
    _previousTranslation = _translation;
    _previousRotation = _rotation;

    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);
    int index;
    float trans[3];
    int previousIndex;
    fscanf(file,"%i %f %f %f %i\n", &index, &trans[0], &trans[1], &trans[2], &previousIndex);

    QUATERNION rotation;
    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
    rotation.x() = quat[0]; rotation.y() = quat[1]; 
    rotation.z() = quat[2]; rotation.w() = quat[3];
    rotation.normalize();
    float length, radius;
    fscanf(file,"%f %f\n", &length, &radius);

    _translation = VEC3F(trans[0], trans[1], trans[2]) * scale;
    _rotation = rotation;

    QUATERNION relativeRotation = _rotation * _restRotation.conjugate();

    VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_restTranslation);

    _dqFromRest = DUAL_QUATERNION(relativeRotation, relativeTranslation);
  }

  void updateTransformAccordingToParent()
  {
    if(_parent == NULL)
      return;
    VEC3F diff = _restTranslation - _parent->restTranslation();
    VEC3F rotatedDiff = _parent->rotation() * _parent->restRotationInv() * diff;
    _translation = _parent->translation() + rotatedDiff;

    QUATERNION relativeRotation = _rotation * _restRotation.conjugate();

    VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_restTranslation);

    _dqFromRest = DUAL_QUATERNION(relativeRotation, relativeTranslation);
  }

  void updateAccordingToParent()
  {
    cout << __FILE__ << " " << __FUNCTION__ << " unimplemented " << endl;
  }

  void interpolateFromRestToRecord(Real step, int recordID)
  {
    if(_parent == NULL){
      _translation = (_restTranslation * (1 - step)) + (_recordedTranslations[recordID] * step);
      QUATERNION fromRot(_restRotation);
      QUATERNION toRot(_recordedRotations[recordID]);
      QUATERNION inbetweenRot = fromRot.slerp(step, toRot);
      inbetweenRot.normalize();
      _rotation = inbetweenRot;
    }else{
      QUATERNION myRelativeRotation(_restRotation * (_parent->_restRotation).conjugate());
      QUATERNION otherRelativeRotation(_recordedRotations[recordID] * (_parent->_recordedRotations[recordID]).conjugate());
      QUATERNION newRelativeRotation = myRelativeRotation.slerp(step, otherRelativeRotation);

      _rotation = newRelativeRotation * (_parent->_rotation).conjugate();
      updateTransformAccordingToParent();
      // VEC3F myRelativeTranslation = _restRotation.transpose() * (_restTranslation - _parent->restTranslation());
    }
  }

  void updateRecentRecordTransformAccordingToParent()
  {
    if(_parent == NULL)
      return;
    VEC3F diff = _restTranslation - _parent->restTranslation();
    VEC3F rotatedDiff = (_parent->_recordedRotations.back() * (_parent->_restRotation).conjugate())._transformVector(diff);
    _recordedTranslations.back() = _parent->_recordedTranslations.back() + rotatedDiff;
  }
  
  inline void setParentID(int id) { _parentId = id; };
  inline void setParent(ODE_BONE* bone) { _parent = bone; }; 
  
  inline VEC3F& translation() { return _translation; };
  inline MATRIX3 rotation() { return _rotation.toRotationMatrix(); };
  // MATRIX4& transform() { return _transform; };

  inline VEC3F& restTranslation() { return _restTranslation; };
  inline MATRIX3 restRotation() { return _restRotation.toRotationMatrix(); };
  inline MATRIX3& restRotationInv() { return _restRotationInv; };

  void reset() {
    _translation = _restTranslation;
    _rotation = _restRotation;
    // _transform = _restTransform;
  };

  inline pair<VEC3F, VEC3F> restBoneSegments()
  {
    VEC3F beginVertex = restRotation() * VEC3F(0, 0, _length * 0.5) + _restTranslation;
    VEC3F endVertex = restRotation() * VEC3F(0, 0, -_length * 0.5) + _restTranslation;
    return make_pair(beginVertex, endVertex);
  }

  void draw() {
    VEC3F beginVertex = rotation() * VEC3F(0, 0, _length * 0.5) + _translation;
    VEC3F endVertex = rotation() * VEC3F(0, 0, -_length * 0.5) + _translation;
    VEC3F centerVertex = _translation;

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
    VEC3F ret = _restRotationInv * (pos - _restTranslation);
    ret = rotation() * ret + _translation;
    return ret;
  }
  inline VEC3F transformFromPrevious(const VEC3F& pos){
    VEC3F ret = _previousRotation.toRotationMatrix().transpose() * (pos - _previousTranslation);
    ret = rotation() * ret + _translation;
    return ret;
  }

  inline void transformFromRestUsingRecords(const VEC3F& pos, vector<VEC3F>& output)
  {
    output.clear();
    VEC3F pullBack = _restRotationInv * (pos - _restTranslation);
    for(unsigned int x = 0; x < _recordedTranslations.size(); x++){
      output.push_back(_recordedRotations[x].toRotationMatrix() * pullBack + _recordedTranslations[x]);
    }
  }
  
  inline DUAL_QUATERNION dualQuaternionRecordFromRest(int recordID)
  {
    assert(recordID >= 0 && recordID < _recordedRotations.size());
    QUATERNION relativeRotation = _recordedRotations[recordID] * _restRotation.conjugate();

    VEC3F relativeTranslation = _recordedTranslations[recordID] - relativeRotation._transformVector(_restTranslation);
    return DUAL_QUATERNION(relativeRotation, relativeTranslation);
  }

  inline DUAL_QUATERNION dualQuaternionFromRest()
  {
    // QUATERNION relativeRotation = _rotation * _restRotation.conjugate();

    // VEC3F relativeTranslation = _translation - relativeRotation._transformVector(_restTranslation);
    // return DUAL_QUATERNION(relativeRotation, relativeTranslation);
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
    // _previousTranslation = _translation;
    // _previousRotation = _rotation;

    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);
    int index;
    float trans[3];
    int previousIndex;
    fscanf(file,"%i %f %f %f %i\n", &index, &trans[0], &trans[1], &trans[2], &previousIndex);

    QUATERNION rotation;
    float quat[4];
    fscanf(file,"%f %f %f %f\n", &quat[0], &quat[1], &quat[2], &quat[3]);
    rotation.x() = quat[0]; rotation.y() = quat[1]; 
    rotation.z() = quat[2]; rotation.w() = quat[3];
    rotation.normalize();
    float length, radius;
    fscanf(file,"%f %f\n", &length, &radius);

    _recordedTranslations.push_back(VEC3F(trans[0], trans[1], trans[2]) * scale);
    _recordedRotations.push_back(rotation);
  };

  inline void clearRecords(){
    _recordedTranslations.clear();
    _recordedRotations.clear();
  };

  bool isRestPose(){
    QUATERNION diffRot = _rotation - _restRotation;
    return (_translation - _restTranslation).norm() + diffRot.x() * diffRot.x() + diffRot.y() * diffRot.y() + diffRot.z() * diffRot.z() < 1e-5;
  }

  bool isRecordPose(int recordID){
    if(recordID < 0 && recordID >= _recordedRotations.size())
      return false;
    QUATERNION diffRot = _rotation - _recordedRotations[recordID];
    return (_translation - _recordedTranslations[recordID]).norm() + diffRot.x() * diffRot.x() + diffRot.y() * diffRot.y() + diffRot.z() * diffRot.z() < 1e-5;
  }

  void computeRelativeRT(ODE_BONE* otherBone, VEC3F& relativeTranslation, QUATERNION& relativeRotation)
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
  VEC3F _translation;
  QUATERNION _rotation;
  VEC3F _previousTranslation;
  QUATERNION _previousRotation;
  // MATRIX4 _transform;

  DUAL_QUATERNION _dqFromRest;

  vector<VEC3F>       _recordedTranslations;
  vector<QUATERNION>  _recordedRotations;

  VEC3F _restTranslation;
  QUATERNION _restRotation;

  MATRIX3 _restRotationInv;

  Real _length;
  Real _radius;
  // MATRIX4 _restTransform;
  ODE_BONE* _parent;
  int _id;
  int _parentId;
};

#endif
