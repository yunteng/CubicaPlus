/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 *****************************************************************************/

#ifndef IO_H
#define IO_H

#include <SETTINGS.h>

#include <dirent.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#if defined(__unix__) || defined (__LINUX__)
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif

using namespace std;

#define SDUMP(x)  " " << #x << "=[ " << x << " ] "
class IO{
public:
  static string intToString(int n)
  {
    stringstream ss;
    ss << n;
    return ss.str().c_str();
  }
  static void split(const string &s, char delim, vector<string>& elems)
  {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
      elems.push_back(item);
    }
  }
  static int getdir(string dir, vector<string> &files)
  {
    DIR *dp;
    struct dirent *dirp;
    if((dp = opendir(dir.c_str())) == NULL) {
      return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
      files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
  }
  ////////////////////////////////////////////////////////////////
  // Print integer to a zero-padded string
  //////////////////////////////////////////////////////////////////
  static std::string itoPaddedString(int frame)
  {
    char buffer[256];
    sprintf(buffer, "%i", frame);

    std::string number = std::string(buffer);
    if (frame < 10) number = std::string("0") + number;
    if (frame < 100) number = std::string("0") + number;
    if (frame < 1000) number = std::string("0") + number;

    return number;
  }
  static vector<string> getAllSnapshots(const string& dataPath){
    cout << "+++++++++++++++++++++++++++++++++" << endl
         << __FILE__ << " " << __FUNCTION__ << endl
         << " Be careful this function only looks for *.state files" << endl
         << "+++++++++++++++++++++++++++++++++" << endl;
    vector<string> filenames;
    vector<string> snapshotIdx;
    if(IO::getdir(dataPath, filenames) < 0){
      cout << "can't get files in directory " << dataPath << "!!!" << endl;
      return snapshotIdx;
    }
      
    for(unsigned int x = 0; x < filenames.size(); x++){
      vector<string> elems;
      IO::split(filenames[x], '.', elems);
      if(elems.size() == 2 && elems[1].compare("state") == 0)
        snapshotIdx.push_back(elems[0]);
    }
    sort(snapshotIdx.begin(), snapshotIdx.end());
    cout << "found " << snapshotIdx.size() << " snapshots" << endl;
    return snapshotIdx;
  }
  // FIXME::very inefficient
  static void write(const VECTOR& vec, const string& filename)
  {
    FILE* file = fopen(filename.c_str(), "wb");
    if (file == NULL)
    {
      cout << " Could not open file " << filename << " to write!" << endl;
      return;
    }
    write(vec, file);
    fclose(file);
  }
  static void write(const VECTOR& vec, FILE* file)
  {
    int size = vec.size();
    fwrite((void*)&size, sizeof(int), 1, file);
    for(int x = 0; x < size; x++){
      double val = vec[x];
      fwrite((void*)&val, sizeof(double), 1, file);
    }
  }
  // FIXME::very inefficient
  static bool read(VECTOR& vec, const string& filename)
  {
    FILE* file;
    file = fopen(filename.c_str(), "rb");
    if (file == NULL)
    {
      cout << " Could not open file " << filename << " to read!" << endl;
      return false;
    }
    int size = 0;
    // read dimensions
    fread((void*)&size, sizeof(int), 1, file);

    if(vec.size() != size)
      vec.resize(size);

    for(int x = 0; x < size; x++){
      double val = 0;
      fread((void*)&val, sizeof(double), 1, file);
      vec[x] = val;
    }
    fclose(file);
    return true;
  }
  static void write(const MATRIX& mat, const string& filename)
  {
    FILE* file = fopen(filename.c_str(), "wb");
    if (file == NULL)
    {
      cout << " Could not open file " << filename << " to read!" << endl;
    }
    write(mat, file);

    fclose(file);
  }
  static void write(const MATRIX& mat, FILE* file)
  {
    int rows = mat.rows();
    int cols = mat.cols();
    double* data = new double[rows * cols];
    int idx = 0;
    for(int x = 0; x < rows; x++)
      for(int y = 0; y < cols; y++)
        data[idx++] = mat(x, y);

    fwrite((void*)&rows, sizeof(int), 1, file);
    fwrite((void*)&cols, sizeof(int), 1, file);
    fwrite((void*)data, sizeof(double), rows * cols, file);

    delete[] data;
  }
  static bool read(MATRIX& mat, const string& filename)
  {
    FILE* file = fopen(filename.c_str(), "rb");
    if(file == NULL)
    {
      cout << "could not open file " << filename << " to read!" << endl;
      return false;
    }
    read(mat, file);
    fclose(file);
    return true;
  }

  static void read(MATRIX& mat, FILE* file)
  {
    int rows = 0;
    int cols = 0;
    fread((void*)&rows, sizeof(int), 1, file);
    fread((void*)&cols, sizeof(int), 1, file);
    double* data = new double[rows * cols];

    fread((void*)data, sizeof(double), rows * cols, file);

    mat.resize(rows, cols);

    int idx = 0;
    for(int x = 0; x < rows; x++)
      for(int y = 0; y < cols; y++)
        mat(x, y) = data[idx++];

    delete[] data;
  }
};

#endif
