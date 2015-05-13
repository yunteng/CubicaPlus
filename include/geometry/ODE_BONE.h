#ifndef ODE_BONE_H
#define ODE_BONE_H

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <SETTINGS.h>
#include <vector>
#include <util/SIMPLE_PARSER.h>
#include <geometry/DUAL_QUATERNION.h>

using namespace::std;

class ODE_BONE
{
public:
  
  ODE_BONE(FILE* file):
    _parent(NULL),
    _parentId(-1)
  {
    Real scale = SIMPLE_PARSER::getFloat("skeleton scaling", 1.0);

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
  
  inline void setParentID(int id) { _parentId = id; };
  inline void setParent(ODE_BONE* bone) { _parent = bone; }; 
  
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

  inline DUAL_QUATERNION& dualQuaternionFromRest()
  {
    return _dqFromRest;
  }

  /*
  compute the relative rotation and translation between this and otherBone
  */
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

  DUAL_QUATERNION _dqFromRest;

  VEC3F _restTranslation;
  QUATERNION _restRotation;

  MATRIX3 _restRotationInv;

  Real _length;
  Real _radius;

  ODE_BONE* _parent;
  int _id;
  int _parentId;
};

#endif
