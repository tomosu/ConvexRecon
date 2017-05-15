#include "cnpy.h"
#include "libarg.h"
#include "StringUtil.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <dirent.h>
#include <cassert>

class BinImageDirectoryOpener{
public:
  BinImageDirectoryOpener(std::string dir_in)
    :dir(dir_in)
  {
    allInputFilePath =getSubDirPath(dir);
    allInputFileNames =getSubDirFilenames(dir);
  }

  int getFileNum(){
    return allInputFilePath.size();
  }

  std::string getFilePath(int index){
    return allInputFilePath[index];
  }

  std::string getFileName(int index){
    return allInputFileNames[index];
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

  std::vector<std::string> getSubDirFilenames(std::string dirPath){
    std::string sla("/");
    std::vector<std::string> subDirFilenames;
    DIR* dp =opendir(dirPath.c_str());
    struct dirent* dent;
    while((dent = readdir(dp)) != NULL){
      std::string filename(dent->d_name);
      if(filename != std::string(".") &&
	 filename != std::string("..") &&
	 filename != std::string(".DS_Store")){
	subDirFilenames.push_back(filename);
      }
    }
    closedir(dp);
    return std::move(subDirFilenames);
  }

  std::string dir;
  std::vector<std::string> allInputFilePath;
  std::vector<std::string> allInputFileNames;
};


std::string extractPrefix(std::string filename){
  std::string dot =".";
  std::vector<std::string> tokens =my_split(filename, dot);
  std::string ret =tokens[0];
  return ret;
}


void convertBin2Npy(std::string input,
		    unsigned int nx, unsigned int ny,
		    std::string output){

  const unsigned int shape[] = {ny,nx};
  std::vector<float> in(nx *ny);
  FILE *fp =NULL;
  if(fp =fopen(input.c_str(), "rb")){
    fread(&in[0], sizeof(float), nx*ny, fp);
    fclose(fp);
  }else{
    std::cout << "invalid input file: " << input << std::endl;
  }
  cnpy::npy_save(output.c_str(), &in[0], shape, 2, "w");
}


#define MBF 256
int main(int argc, char *argv[])
{
  std::string sla("/");
  std::string npy(".npy");
  int nx =512;
  int ny =512;
  char inputDir[MBF] ="./output";
  char outputDir[MBF] ="./npy_output";

  ARGUMENT arg;
  arg.addOption("-input_dir", inputDir, MBF);
  arg.addOption("-output_dir", outputDir, MBF);

  if(0 != arg.parse(argc, argv) || argc !=1){
    printf("wrong option parameters\n");
    exit(1);
  }

  BinImageDirectoryOpener inputFiles( inputDir );

  for(int i=0; i<inputFiles.getFileNum(); i++){
    std::string outFilepath = std::string(outputDir) +sla
      +extractPrefix(inputFiles.getFileName(i)) +npy;
    convertBin2Npy(inputFiles.getFilePath(i),
		   nx, ny, outFilepath);
  }
}
