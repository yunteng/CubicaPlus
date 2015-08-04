/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef DUAL_QUATERNION_H
#define DUAL_QUATERNION_H

#include <SETTINGS.h>
#include <util/MATRIX_UTIL.h>

class DUAL_QUATERNION
{
public:
  DUAL_QUATERNION()
  {
    // eigen's quaternion's constructor is (w, x, y, z)
    // storage order is x, y, z, w
    real = QUATERNION(1, 0, 0, 0);
    dual = QUATERNION(0, 0, 0, 0);
  }
  DUAL_QUATERNION(const QUATERNION& r, const QUATERNION& d)
  {
    real = r;
    dual = d;
  }
  DUAL_QUATERNION(const QUATERNION& q, const VEC3F& t)
  {
    real = q;
    dual = (QUATERNION(0, t[0], t[1], t[2]) * real) * 0.5;
  }
  inline Real dot(const DUAL_QUATERNION& other) const
  {
    return real.dot(other.real);
  }
  inline DUAL_QUATERNION operator*(Real scale) const
  {
    return DUAL_QUATERNION(real * scale, dual * scale);
  }
  inline DUAL_QUATERNION& normalize()
  {
    Real mag = real.norm();
    assert(mag > 1e-6);
    real = real * (1.0 / mag);
    dual = dual * (1.0 / mag);
    return *this;
  }
  inline DUAL_QUATERNION normalized() const
  {
    Real mag = real.norm();
    assert(mag > 0);
    return DUAL_QUATERNION(real * (1.0 / mag), dual * (1.0 / mag));
  }

  inline DUAL_QUATERNION operator+(const DUAL_QUATERNION& other) const
  {
    return DUAL_QUATERNION(real + other.real, dual + other.dual);
  }

  // Multiplication order - left to right
  inline DUAL_QUATERNION operator*(const DUAL_QUATERNION& other) const
  {
    return DUAL_QUATERNION(other.real * real, (other.dual * real) + (other.real * dual));
  }

  inline DUAL_QUATERNION conjugate()
  {
    return DUAL_QUATERNION(real.conjugate(), dual.conjugate());
  }
  
  inline QUATERNION rotation() const
  {
    return real;
  }
  inline VEC3F translation() const
  {
    QUATERNION t = (dual * 2.0) * real.conjugate();
    return VEC3F(t.x(), t.y(), t.z());
  }
  inline VEC3F transform(const VEC3F& vec) const
  {
    return real._transformVector(vec) + translation();
  }
  inline VEC3F inverseTransform(const VEC3F& vec) const
  {
    return real.conjugate()._transformVector(vec - translation());
  }

public:
  QUATERNION real;
  QUATERNION dual;


};

#endif
