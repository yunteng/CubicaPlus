/*
This file is part of Cubica.
 
Cubica is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Cubica is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Cubica.  If not, see <http://www.gnu.org/licenses/>.
*/
// 
// simplistic config file parser
//
//////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cctype>
#include <util/SIMPLE_PARSER.h>
#include <material/COROTATION.h>

//////////////////////////////////////////////////////////////////////////////
// whitespace helper
//////////////////////////////////////////////////////////////////////////////
static string removeWhitespace(string in) {
  size_t s = 0;
  size_t e = in.length();
  while(s<in.length() && isspace(in[s])) { s++; }
  while(e>s && isspace(in[e-1])) { e--; }
  return in.substr(s,e-s);
}

//////////////////////////////////////////////////////////////////////////////
// force string to lowercase
//////////////////////////////////////////////////////////////////////////////
static void forceLower(string& input) {
  string::iterator i;
  for (i = input.begin(); i != input.end(); i++)
    *i = tolower(*i);
}
map<string, string> SIMPLE_PARSER::mVals;
map<string, bool> SIMPLE_PARSER::mUsed;

//////////////////////////////////////////////////////////////////////////////
// Constructor / Destructor
//////////////////////////////////////////////////////////////////////////////
bool SIMPLE_PARSER::parse(std::string file)
{
  if(file.length()==0) {
    std::cout << "Skipping config file\n";
    return false;
  }

  std::cout << "Using config file "<< file <<"\n";

  int lineCnt=1;

  char buffer[512];
  std::ifstream myfile (file.c_str());

  if (!myfile.good())
  {
    cout << " Failed to open file " << file << "!!!" << endl;
    return false;
  }
  if (myfile.is_open())
  {
    while (! myfile.eof() )
    {
      myfile.getline(buffer, 512);
      string line(buffer);
      if(line.length()<1) continue;
      if(line[0] == '#') continue;

      size_t pos = line.find_first_of("=");
      if(pos != string::npos && pos<line.length() ) {
        string lhs = removeWhitespace( line.substr(0,pos) );
        string rhs = removeWhitespace( line.substr(pos+1,line.length()-pos) );

        forceLower(lhs);
        forceLower(rhs);

        // store...
        mVals[lhs] = rhs;
      } else {
        // simple check for errors...
        string check = removeWhitespace(line);
        if(check.length()>0) {
          std::cerr<<"Unable to parse, error in line "<<lineCnt<<": '"<< line <<"' !\n";
          // exit(1);
          return false;
        }
      }
      lineCnt++;
    }
    myfile.close();
    return true;
  } 
  else 
  {
    std::cerr<<"Unable to parse!\n";
    return false;
  }
}


//////////////////////////////////////////////////////////////////////////////
// See if a parameter was defined
//////////////////////////////////////////////////////////////////////////////
bool SIMPLE_PARSER::defined(string name)
{
  map<string, string>::iterator i;
  i = mVals.find(name);

  return (i != mVals.end());
}

//////////////////////////////////////////////////////////////////////////////
// generic scalar retrieval
//////////////////////////////////////////////////////////////////////////////
template<class T> T SIMPLE_PARSER::getScalarValue(string name, T defaultValue, bool needed)
{
  T ret = 0;
  forceLower(name);
  if(mVals.find(name) == mVals.end()) {
    if(needed) {
      std::cerr<<"Required value '"<<name<<"' not found in config file!\n";
      exit(1); 
    }
    return defaultValue;
  }
  ret = (T)atof(mVals[name].c_str());
  mUsed[name] = true;
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// get an VEC3F 
//////////////////////////////////////////////////////////////////////////////
VEC3F SIMPLE_PARSER::getVEC3F(string name, VEC3F defaultValue, bool needed)
{
  VEC3F ret;
  forceLower(name);
  if(mVals.find(name) == mVals.end()) {
    if(needed) {
      std::cerr << "Required value '" << name << "' not found in config file!\n";
      exit(1); 
    }
    return defaultValue;
  }
  vector<string> tokens;
  IO::split(mVals[name], ',', tokens);
  for(int v = 0; v < 3; v++)
    ret[v] = atof(tokens[v].c_str());
  mUsed[name] = true;
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// get an integer 
//////////////////////////////////////////////////////////////////////////////
int SIMPLE_PARSER::getInt(string name,    int defaultValue,    bool needed)
{
  int ret = getScalarValue<int>(name, defaultValue, needed);
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// get a boolean
//////////////////////////////////////////////////////////////////////////////
bool SIMPLE_PARSER::getBool(string name, bool defaultValue, bool needed)
{
  bool ret = (getScalarValue<int>(name, defaultValue, needed) != 0);
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// get a floating point
//////////////////////////////////////////////////////////////////////////////
double SIMPLE_PARSER::getFloat(string name,  double defaultValue, bool needed)
{
  double ret = getScalarValue<double>(name, defaultValue, needed);
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// get a string
//////////////////////////////////////////////////////////////////////////////
string SIMPLE_PARSER::getString(string name, string defaultValue, bool needed)
{
  string ret("");
  forceLower(name);
  if(mVals.find(name) == mVals.end()) {
    if(needed) {
      std::cerr<<"Required value '"<<name<<"' not found in config file!\n";
      exit(1); 
    }
    return defaultValue;
  }
  ret = mVals[name];
  mUsed[name] = true;

  // force to lower case
  forceLower(ret);

  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// check if there were any unused pairs
//////////////////////////////////////////////////////////////////////////////
bool SIMPLE_PARSER::haveUnusedValues()
{
  for(std::map<string, string>::iterator i=mVals.begin();
      i != mVals.end(); i++) {
    if((*i).second.length()>0) {
      if(!mUsed[ (*i).first]) {
        return true;
      }
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// print unused pairs
//////////////////////////////////////////////////////////////////////////////
string SIMPLE_PARSER::printAllUnused()
{
  std::ostringstream out;
  for(std::map<string, string>::iterator i=mVals.begin();
      i != mVals.end(); i++) {
    if((*i).second.length()>0) {
      if(!mUsed[ (*i).first ]) {
        out <<"'"<<(*i).first<<"'='"<<(*i).second<<"' ";
      }
    }
  }
  return out.str();
}

//////////////////////////////////////////////////////////////////////////////
// Read a material file and return a pointer to it
//////////////////////////////////////////////////////////////////////////////
MATERIAL *SIMPLE_PARSER::readMaterial(string filename)
{
  parse(filename);
  // set material
  MATERIAL* material = NULL;
  string materialType("");
  materialType = getString("material type", materialType);

  /*if (materialType.compare("stvk") == 0)
  {
    double lambda = 10.0;
    double mu = 50.0;
    lambda = getFloat("stvk lambda", lambda);
    mu = getFloat("stvk mu", mu);
    material = new STVK(lambda, mu);
    cout << "==================================================" << endl;
    cout << " Material is St-VK, lambda = " << lambda << " mu = " << mu << endl;
  }else*/
  if (materialType.compare("corotation") == 0)
  {
    double lambda = 10.0;
    double mu = 50.0;
    bool computeLame = false;
    computeLame = getBool("compute lame", computeLame);
    if(!computeLame){
      lambda = getFloat("corotation lambda", lambda);
      mu = getFloat("corotation mu", mu);
    }else{
      double young = getFloat("corotation young", 4000);
      double possion = getFloat("corotation possion", 20000);
      lambda =  young * possion / ( 1.0 + possion) / ( 1.0 - 2.0 * possion);
      mu = young / 2.0 / (1.0 + possion);
      // cout << "lambda " << lambda << " mu " << mu << endl;
    }
    cout << "==================================================" << endl;
    material = new COROTATION(lambda, mu);
    cout << " Material is corotation, lambda = " << lambda << " mu = " << mu << endl;
    cout << "==================================================" << endl;
  }
  /*else if (materialType.compare("mooney-rivlin") == 0)
  {
    double mu01 = 100.0;
    double mu10 = 500.0;
    double k = 100000.0;
    mu01 = getFloat("mooney-rivlin mu01", mu01);
    mu10 = getFloat("mooney-rivlin mu10", mu10);
    k = getFloat("mooney-rivlin k", k);
    material = new MOONEY_RIVLIN(mu01, mu10, k);
    cout << "==================================================" << endl;
    cout << " Material is Mooney-Rivlin, mu01 = " << mu01 << " mu10 = " << mu10 << " k = " << k << endl;
  }
  else if (materialType.compare("arruda-boyce") == 0)
  {
    double nkTheta = 5000.0;
    double N  = 5.0;
    double k = 1000.0;
    nkTheta = getFloat("arruda-boyce nktheta", nkTheta);
    N = getFloat("arruda-boyce n", N);
    k = getFloat("arruda-boyce k", k);
    material = new ARRUDA_BOYCE(nkTheta, N, k);
    cout << "==================================================" << endl;
    cout << " Material is Arruda-Boyce, nkTheta = " << nkTheta << " N = " << N << " k = " << k << endl;
  }
  else if (materialType.compare("neo-hookean") == 0)
  {
    double mu = 50.0;
    double lambda = 10.0;
    mu = getFloat("neo-hookean mu", mu);
    lambda = getFloat("neo-hookean lambda", lambda);
    material = new NEO_HOOKEAN(mu, lambda);
    cout << "==================================================" << endl;
    cout << " Material is Neo-Hookean, mu = " << mu << " lambda = " << lambda << endl;
  }*/
  else
  {
    cout << " *** Material type undefined! *** " << endl;
    exit(1);
  }

  // set the invertible wrapper if necessary
  // bool invertible = false;
  // invertible = SIMPLE_PARSER::getBool("invertible", invertible);
  // if (invertible)
  // {
  //   cout << " Setting material to invertible" << endl;
  //   material = new INVERTIBLE(material);
  // }
  return material;
}

void SIMPLE_PARSER::setString(string name, string val)
{
  forceLower(name);
  mVals[name] = val;
}

