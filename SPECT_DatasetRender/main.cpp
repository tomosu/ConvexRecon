/************************************************************
 *
 * Basic 2D-MLEM solver
 *
 ************************************************************/
#include "libini.h"
#include "libarg.h"

#include "Geometry.hpp"
#include "Generator.h"
#include "FP.hpp"
#include "FFT.h"
#include "Phantom.h"
#include "GeometricTools.h"
#include "ProjectionGenerator.h"
#include "ProjectionOperator.hpp"
#include "DicomSliceReader.hpp"
#include "DicomDirectryOpener.hpp"
#include "StringUtil.h"
#include "NpyReconConverter.hpp"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <algorithm>
#include <complex>
#include <random>

#define MBF 256



void AddCountNoise(const Geometry &geo,
		   Projection  *proj,
		   unsigned int seed){

  std::mt19937 mt(seed);
  std::uniform_real_distribution<double> genrand1(0.0, 1.0);

  for (int prjID = 0; prjID < geo.projectionCount; prjID++) {
    for (int detID = 0; detID < geo.detectorPixelNum; detID++) {
      while (1) {
	// BoxMuller
	Real_t val =proj->getValue(detID, prjID);
	Real_t SND   = sqrtf(-2 * logf(genrand1(mt)) ) * cosf(2 * M_PI * genrand1(mt));
	Real_t count = val + sqrt(val) * SND;
	if (count >= 0) {
	  proj->setValue(detID, prjID,count);
	  break;
	} else {
	  continue;
	}
      }
    }
  }

}


/****************************************
 * Shepp-Logan FBP
 */
class DirectryOpener{

public:
  DirectryOpener(std::string rootDir_in, int level)
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

    for(auto f: allFilePath){
      allFileNames.push_back(CreateOmittedFilename(f));
    }
  }

  int getFileNum(){
    return allFilePath.size();
  }

  std::string getFilePath(int index){
    return allFilePath[index];
  }

  std::string getFileName(int index){
    return allFileNames[index];
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


  std::string CreateOmittedFilename(std::string filepath){
    std::string sla ="/";
    std::vector<std::string> tokens =my_split(filepath, sla);
    return tokens[tokens.size()-1];
  }

  std::string rootDir;
  std::string prefix;
  std::vector<std::string> allFilePath;
  std::vector<std::string> allFileNames;
};






/****************************************
 * Shepp-Logan FBP
 */
inline void FT_1d(int ir,
		  std::vector< std::complex<Real_t> > in,
		  std::vector< std::complex<Real_t> > &out)
{
  int nx =in.size();
  Real_t  t  = 2 *M_PI /nx;
  Real_t *gr = (Real_t *)malloc((unsigned long)nx * sizeof(Real_t));
  Real_t *gi = (Real_t *)malloc((unsigned long)nx * sizeof(Real_t));

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < nx; i++) {
    Real_t u = i -nx /2;
    gr[i] = gi[i] = 0;
    for (int j = 0; j < nx; j++) {
      Real_t x = j - nx / 2;
      gr[i] += (Real_t)(in[j].real() *cos(t *u *x)        +ir *in[j].imag() *sin(t *u *x));
      gi[i] += (Real_t)(in[j].real() *sin(t *u *x) *(-ir) +in[j].imag() *cos(t *u *x));
    }
  }

  int n =1;
  if (ir == -1) n =nx;

  for (int i =0; i <nx; i++) {
    out[i].real(gr[i] /n);
    out[i].imag(gi[i] /n);
  }

  free(gr);
  free(gi);
}


Reconstruction FBP(const Geometry &geo,
		   const Projection &data)
{
  ProjectionOperator op(geo);

  // set filter
  std::vector< std::complex<Real_t> > shepp(geo.detectorPixelNum);
  {
    for (int i = 0; i < geo.detectorPixelNum; i++) {
      int n = abs(i - geo.detectorPixelNum / 2);
      shepp[i] =std::complex<Real_t>(2.0 /M_PI /M_PI /(1.0 -4.0 *n *n), 0.0);
    }
  }
  std::vector< std::complex<Real_t> > shepp_f(geo.detectorPixelNum);
  FT_1d(1, shepp, shepp_f);

  // set filtered projection
  Projection filteredProjection =data;
  std::cout << "f init proj stat:" << std::endl;
  filteredProjection.stat();
  for (int i = 0; i < geo.projectionCount; i++) {

    std::vector< std::complex<Real_t> > view(geo.detectorPixelNum);
    for (int j = 0; j < geo.detectorPixelNum; j++) {
      view[j].real( data.getValue(j,i) );
      view[j].imag( 0.0 );
    }

    std::vector< std::complex<Real_t> > view_f(geo.detectorPixelNum);
    FT_1d(1, view, view_f);
    for (int j = 0; j < geo.detectorPixelNum; j++) {
      view_f[j] = view_f[j] *shepp_f[j];
    }
    FT_1d(-1, view_f, view);

    for (int j = 0; j < geo.detectorPixelNum; j++) {
      filteredProjection.setValue(j, i, view[j].real());
    }
  }
  std::cout << "f proj stat:" << std::endl;
  filteredProjection.stat();

  //BP
  Reconstruction ret =op.backwardProjection(filteredProjection);
  return std::move(ret);
}



/****************************************
 * util name
 */
std::string CreateOmittedFilename(std::string filepath){
  std::string sla ="/";
  std::vector<std::string> tokens =my_split(filepath, sla);
  std::string ret =tokens[4] +std::string("_") +tokens[7];
  return ret;
}



/****************************************
 * main
 */
int main(int argc, char **argv)
{

  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;

  // arg
  ARGUMENT arg;
  #define MBF 256
  char outputDir[MBF] ="output";
  double emitRatio =0.1;
  double attenRatio =0.0001;
  arg.addOption("-output_dir", outputDir, MBF);
  arg.addOption("-emit_ratio", emitRatio);
  arg.addOption("-atten_ratio", attenRatio);
  int retcd = arg.parse(argc, argv);
  if (retcd || argc != 1) {
    printf("Usage: %s\n", argv[0]);
    return -1;
  }



  // general
  char   geometryPath[MBF]           = "../config/geometry_128.ini";
  double noiseIntensity              = 0.25;
  Geometry geo(geometryPath);

  auto petFiles =DirectryOpener(std::string("../Registration/reg_images/pet"), 1);
  auto ctFiles =DirectryOpener(std::string("../Registration/reg_images/ct"), 1);

  int studyNum =petFiles.getFileNum();
  for(int i=0; i<studyNum; i++){

    std::cout << "current image index:" << i << std::endl;
    /****************************************
     * set answer
     */
    std::string filename =petFiles.getFileName(i);
    Reconstruction emitIm =Npy2Reconstruction(petFiles.getFilePath(i).c_str(), geo);
    Reconstruction attenIm =Npy2Reconstruction(ctFiles.getFilePath(i).c_str(), geo);
    attenIm =attenIm.add_scalar(1000);
    std::cout << "emit:" << std::endl;
    emitIm.stat();
    std::cout << "attenuation:" << std::endl;
    attenIm.stat();

    /****************************************
     * generate projection
     */
    //Real_t  emitRatio =0.1;
    //Real_t  emitRatio =0.1;
    //Real_t attenRatio =0.0001;
    //Real_t attenRatio =0.0006;

    EmissionProjectionOperator epo(geo);
    Projection emitProj =epo.forwardProjectionWithAttenuation(emitIm, attenIm,
							      emitRatio, attenRatio);
    int seed =0;
    AddCountNoise(geo, &emitProj, seed);
    std::cout << "emit proj:" << std::endl;
    emitProj.stat();


    /****************************************
     * FBP
     */
    {
      char outPath[MBF];
      sprintf(outPath, "./%s/%s", outputDir, filename.c_str());
      Reconstruction fbpRec =FBP(geo, emitProj);
      std::cout << "FBP rec stat:" << std::endl;
      fbpRec.stat();
      Reconstruction2Npy(outPath, fbpRec);
    }
  }
}
