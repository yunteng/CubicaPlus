#ifndef TIMING_BREAKDOWN_H
#define TIMING_BREAKDOWN_H

#include <SETTINGS.h>
#include <stdio.h>
#include <iostream>
#include <map>
#include <util/TIMER.h>

using namespace::std;

class TIMING_BREAKDOWN
{
public:
  inline static void startFrame(){
    _totalTimer.restart();
  }
  inline static void endFrame(){
    _totalTime += _totalTimer.timing();
    _frameTime.push_back(_totalTimer.timing());
    _totalSteps++;
  }
  inline static void tic(){
    if(_ticCnt == 0)
      _currentTimer.restart();
    _ticCnt++;
  };
  inline static void toc(string function){
    if(_ticCnt == 1)
      _timingBreakdown[function] += _currentTimer.timing();
    _ticCnt--;
  }
  static void printTimingBreakdown(){
    // create an inverse map so that it will sort by time
    map<double, string> inverseMap;
    
    for(map<string, double>::iterator forwardIter = _timingBreakdown.begin(); forwardIter != _timingBreakdown.end(); forwardIter++)
      inverseMap[forwardIter->second] = forwardIter->first;

    // print the map out backwards since it sorts from least to greatest
    cout << "==============================================================================" << endl;
    cout << "TIMING BREAKDOWN: " << endl;
    cout << "==============================================================================" << endl;
    
    double totalSeen = 0.0;
    stringstream ss;
    for (map<double,string>::reverse_iterator backwardIter = inverseMap.rbegin(); backwardIter != inverseMap.rend(); backwardIter++)
    {
      string name = (*backwardIter).second;
      name = name.substr(0,45);
      while (name.size() < 45)
        name = name + string(" ");

      ss.str(" ");
      ss << (*backwardIter).first / _totalTime * 100.0 << "%";
      string percent(ss.str());
      while (percent.size() < 12)
        percent = percent + string(" ");

      //cout << "[" << (*backwardIter).first / _totalTime * 100.0 << "%\t]: "
      cout << "[" << percent.c_str() << "]: "
           << name.c_str() << "\t" << (*backwardIter).first / _totalSteps << "s / frame" << endl;
      totalSeen += (*backwardIter).first;
    }
    Real misc = (_totalTime - totalSeen) / _totalTime * 100.0;

    ss.str(" ");
    ss << misc << "%";
    string percent(ss.str());
    while (percent.size() < 12)
      percent = percent + string(" ");
    cout << "[" << percent << "]: " << "Misc. " << endl;
    cout << "==============================================================================" << endl;
    cout << " Current FPS: " << _totalSteps / _totalTime << endl;
    // cout << " Newton steps per second: " << _totalNewtonStepsSeen / _totalTime << endl;
    // cout << " Mean Newton steps: " << (Real)_totalNewtonStepsSeen / (Real)_totalSteps << endl;
    // cout << " Total Newton stalls: " << _newtonStalls << endl;
    cout << "==============================================================================" << endl;
    if (misc < 0.0)
    {
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << " BREAKDOWN ADDS UP TO MORE THAN 100! TIMERS ARE OVERLAPPING! " << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    }
  }
  static void writeFrameTime(const string& filename)
  {
    FILE* file = fopen(filename.c_str(), "w");
    if(file == NULL){
      cout << __FILE__ << " " << __FUNCTION__ << endl;
      cout << " Cannot open " << filename << " to write!!!" << endl;
      return;
    }

    for(unsigned int x = 0; x < _frameTime.size(); x++){
      fprintf(file, "%f\n", _frameTime[x]);
    }
    fclose(file);
  }
  static void clear()
  {
    _timingBreakdown.clear();
    _totalSteps = 0;
  }

private:
  static map<string, double> _timingBreakdown;
  static vector<Real> _frameTime;
  static TIMER _totalTimer;
  static TIMER _currentTimer;
  static Real _totalTime;
  static int _totalSteps;
  static int _ticCnt;

};

// map<string, double> TIMING_BREAKDOWN::_timingBreakdown;
// TIMER TIMING_BREAKDOWN::_totalTimer;
// TIMER TIMING_BREAKDOWN::_currentTimer;

// Real TIMING_BREAKDOWN::_totalTime = 0;
// int TIMING_BREAKDOWN::_totalSteps = 0;
// // to detect if tic is called multiple times before toc
// int TIMING_BREAKDOWN::_ticCnt = 0;

#endif
