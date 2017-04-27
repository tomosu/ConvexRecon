#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <cassert>

class DicomDirectryOpener{

public:
  DicomDirectryOpener(std::string rootDir_in, int level)
    :rootDir(rootDir_in)
  {
    assert(level > 0);
    std::vector<std::string> sub =getSubDirPath(rootDir);
    for(int i=1; i <level; i++){
      std::vector<std::string> subNext;
      for(auto s : sub){
	for(auto sn : getSubDirPath(s)){
	  subNext.push_back(sn);
	}
      }
      sub =subNext;
    }
    allFilePath =sub;
  }

  int getFileNum(){
    return allFilePath.size();
  }

  std::string getFilePath(int index){
    return allFilePath[index];
  }

private:

  std::vector<std::string> getSubDirPath(std::string dirPath){
    std::string sla("/");

    std::vector<std::string> subDirPath;
    DIR* dp =opendir(dirPath.c_str());
    struct dirent* dent;
    while((dent = readdir(dp)) != NULL){
      std::string filename(dent->d_name);
      if(filename != std::string(".") &&
	 filename != std::string("..") &&
	 filename != std::string(".DS_Store")){
	std::string path =dirPath +sla +filename;
	subDirPath.push_back(path);
      }
    }
    closedir(dp);
    return std::move(subDirPath);
  }

  std::string rootDir;
  std::string prefix;
  std::vector<std::string> allFilePath;
  std::vector<std::string> allFileNames;
};
