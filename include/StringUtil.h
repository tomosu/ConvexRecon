#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>




////////////////////////////
// String Util
////////////////////////////
inline double my_stod(const std::string &in)
{
  double ret;
  std::istringstream ss;
  ss.str(in);
  ss >> ret;
  return ret;
}



inline int my_stoi(const std::string &in)
{
  int ret;
  std::istringstream ss;
  ss.str(in);
  ss >> ret;
  return ret;
}




inline std::vector<std::string> my_split(const std::string &str, const std::string &delim)
{
  std::vector<std::string> res;
  size_t current=0, found, delimlen=delim.size();
  while( (found = str.find(delim, current)) != std::string::npos ){
    res.push_back( std::string(str, current, found-current) );
    current = found +delimlen;
  }
  res.push_back( std::string(str, current, str.size()-current) );
  return res;
}

