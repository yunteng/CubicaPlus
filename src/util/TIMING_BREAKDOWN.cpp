/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#include <util/TIMING_BREAKDOWN.h>

map<string, double> TIMING_BREAKDOWN::_timingBreakdown;
vector<Real> TIMING_BREAKDOWN::_frameTime;
TIMER TIMING_BREAKDOWN::_totalTimer;
TIMER TIMING_BREAKDOWN::_currentTimer;

Real TIMING_BREAKDOWN::_totalTime = 0;
int TIMING_BREAKDOWN::_totalSteps = 0;
// to detect if tic is called multiple times before toc
int TIMING_BREAKDOWN::_ticCnt = 0;
