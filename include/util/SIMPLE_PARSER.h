/******************************************************************************
 *
 * Copyright (c) 2015, Yun Teng (yunteng.cs@cs.ucsb.edu), University of
 # California, Santa Barbara.
 * All rights reserved.
 *
 * Original source code was provided courtesy of 
 * Nils Theurey (www.ntoken.com)
 * Modified by Theodore Kim (http://www.mat.ucsb.edu/~kim/) and
 * Yun Teng
 *****************************************************************************/

#ifndef SIMPLE_PARSER_H
#define SIMPLE_PARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <util/IO.h>
#include <material/MATERIAL.h>

using namespace std;

// read the whole test file of pairs
// a = b
// and look up values with conversion
class SIMPLE_PARSER  
{
  public:
    // SIMPLE_PARSER(std::string file);
    static bool parse(std::string file);
    // get values from file
    static int      getInt(string name,    int defaultValue,    bool needed=false);
    static bool     getBool(string name,   bool defaultValue,   bool needed=false);
    static double   getFloat(string name,  double defaultValue, bool needed=false);
    static string   getString(string name, string defaultValue, bool needed=false);
    static VEC3F    getVEC3F(string name,  VEC3F defaultValue,  bool needed=false);
    // check if there were any unused pairs
    static bool haveUnusedValues();
    static string printAllUnused();

    // check if the a parameters was specified in the config file
    static bool defined(string name);

    static MATERIAL* readMaterial(string filename);

    static void setString(string name, string val);
    
  protected:
    // value pair storage
    static map<string, string> mVals;
    // value pair check if used...
    static map<string, bool> mUsed;

    template<class T> T static getScalarValue(string name, T defaultValue, bool needed=false);
    // template<class T> TVEC3<T> static getVEC3Value(string name, TVEC3<T> defaultValue, bool needed=false);
};

#endif

