#ifndef RGB_HSV_H
#define RGB_HSV_H

#include <SETTINGS.h>

class HSV;
class RGB {
public:
  RGB(): r(0), g(0), b(0) {};
  RGB(double red, double green, double blue):
    r(red), g(green), b(blue) {};
  VEC3F toVEC3F() const{ return VEC3F(r, g, b); };
  HSV toHSV() const;
  double r;       // percent
  double g;       // percent
  double b;       // percent
  
};

class HSV {
public:
  HSV(): h(0), s(0), v(0) {};
  HSV(double hue, double sat, double val):
    h(hue), s(sat), v(val) {};
  VEC3F toVEC3F(){ return VEC3F(h, s, v); };
  HSV lerp(const HSV& end, Real t)
  {
    return HSV((1 - t) * h + t * end.h,
               (1 - t) * s + t * end.s,
               (1 - t) * v + t * end.v);
  }
  RGB toRGB() const;
  double h;       // angle in degrees
  double s;       // percent
  double v;       // percent
};

#endif
