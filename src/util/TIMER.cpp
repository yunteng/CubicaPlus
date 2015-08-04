/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 * Original source code was provided courtesy of 
 * Theodore Kim (http://www.mat.ucsb.edu/~kim/)
 * Modified by Yun Teng
 *****************************************************************************/
#include <util/TIMER.h>

#if _WIN32
bool TIMER::_firstCall = true;
double TIMER::_divisor;
#endif
