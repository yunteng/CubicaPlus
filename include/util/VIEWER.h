#ifndef VIEWER_H
#define VIEWER_H

#include <SETTINGS.h>
#include <glvu.hpp>
#include <iostream>
// #include <IO.h>
#include <util/SIMPLE_PARSER.h>
#include <util/TIMING_BREAKDOWN.h>
// #include <QUICKTIME_MOVIE.h>

using namespace std;
template<class T>
class VIEWER
{
public:
  static void init();
  static void screenshot(string renderPath, int frame);
  ~VIEWER() { delete simulator; };

public:
  static T* simulator;
  static GLVU glvu;
  static bool animate;
  static bool step;

private:
  static void displayFunc();
  static void keyboardFunc(unsigned char Key, int x, int y);
  static void idleFunc();
  static void mouseFunc(int button, int state, int x, int y);
  static void motionFunc(int x, int y);
  static int glvuWindow(glvuVec3f bboxCenter, Real eyeDistanceScale = 1.0);
  static VEC3F unproject(float x, float y, float z);
private:
  static const int windowStartX = 0;
  static const int windowStartY = 0;
  static const int windowWidth = 700;
  static const int windowHeight = 700;
};
#include "VIEWER.inl"

#endif
