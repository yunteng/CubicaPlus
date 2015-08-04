/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <geometry/SURFACE.h>

SURFACE::SURFACE():
  _type("abstract surfae"),
  _thickness(0)
{
  _collisionStiffness = SIMPLE_PARSER::getFloat("collision stiffness", 100);
  _collisionDamping = SIMPLE_PARSER::getFloat("collision damping", 0.1);
  cout << " Collision stiffness " << _collisionStiffness << ", collision damping " << _collisionDamping << endl;
}
SURFACE::~SURFACE()
{
  
}
