#include <util/RGB_HSV.h>

HSV RGB::toHSV() const
{
  HSV         out;
  double      min, max, delta;

  min = r < g ? r : g;
  min = min  < b ? min  : b;

  max = r > g ? r : g;
  max = max  > b ? max  : b;

  out.v = max;                                // v
  delta = max - min;
  if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
      out.s = (delta / max);                  // s
  } else {
      // if max is 0, then r = g = b = 0              
          // s = 0, v is undefined
      out.s = 0.0;
      out.h = NAN;                            // its now undefined
      return out;
  }
  if( r >= max )                           // > is bogus, just keeps compilor happy
      out.h = ( g - b ) / delta;        // between yellow & magenta
  else
  if( g >= max )
      out.h = 2.0 + ( b - r ) / delta;  // between cyan & yellow
  else
      out.h = 4.0 + ( r - g ) / delta;  // between magenta & cyan

  out.h *= 60.0;                              // degrees

  if( out.h < 0.0 )
      out.h += 360.0;

  return out;
}

RGB HSV::toRGB() const
  {
    double      hh, p, q, t, ff;
    long        i;
    RGB         out;

    if(s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = v;
        out.g = v;
        out.b = v;
        return out;
    }
    hh = h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s * ff));
    t = v * (1.0 - (s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = v;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = v;
        break;
    case 5:
    default:
        out.r = v;
        out.g = p;
        out.b = q;
        break;
    }
    return out; 
  };
