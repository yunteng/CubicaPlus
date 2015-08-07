/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/
#ifndef NNHTP_CUBATURE_GENERATOR_H
#define NNHTP_CUBATURE_GENERATOR_H

#include <SETTINGS.h>
#include <fstream>
#include <util/MERSENNETWISTER.h>
/*
implementation of the Non-negative Hard-Threshold Pursuing algorithm in C. von Tycowicz, C. Schulz, H.-P. Seidel, and K. Hildebrandt, An efficient construction of reduced deformable objects.
*/
template<class APP>
class NNHTP_CUBATURE_GENERATOR
{
public:
  NNHTP_CUBATURE_GENERATOR(APP* app);
  ~NNHTP_CUBATURE_GENERATOR();

  void generateCubatures();
  void verifyCubaturePoints(const VECTOR& b);
  void writeCubatures(const string& filename);

  vector<int>& keyPointIDs() { return _keyPointIDs; };
  vector<Real>& keyWeights() { return _keyWeights; };

private:
  APP* _app;
  MERSENNETWISTER _trainingTwister;

  vector<int> NNHTP_randomPickCandidates(vector<int>& excludes, int candidatesPerTry, int totalCandidates, bool addExcludes);
  vector<int> NNHTP_project(VECTOR& w, int toKeep);

  vector<int> _keyPointIDs;
  vector<Real> _keyWeights;
  int _nCutoff;

  struct orderByValueSmallerThan{
    bool operator()(pair<int, Real> const& a, pair<int, Real> const& b) const{
      return a.second < b.second;
    }
  };
};

#include "NNHTP_CUBATURE_GENERATOR.inl"
#endif
