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
#include <iostream>
#include <fstream>

#define MBF 256




/****************************************
 * initialize pixel
 */
inline void InitializePixel(const Geometry &geo,
                            const Projection &proj,
			    Reconstruction *recon)
{
  // calculate mean projection value per 1 detector
  Real_t projectionMean = 0.0F;
  for (int i = 0; i < geo.projectionCount; i++) {
    for (int j = 0; j < geo.detectorPixelNum; j++) {
      projectionMean += proj.getValue(j, i);
    }
  }

  // calculate mean atten per 1 pixel
  Real_t initialValue = projectionMean / geo.volumeSize.x; // tekitou
  recon->clear(initialValue);
}



/****************************************
 * Poisson noise on projection
 */
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
 * MLEM Loop
 */
Reconstruction MLEMLoop(const Geometry &geo,
			const Reconstruction &nowRecEmitIm,
			const Reconstruction &attenIm,
			double emitRatio,
			double attenRatio,
			const Projection &data)
{
  EmissionProjectionOperator op(geo);

  //FP
  Projection fproj =op.forwardProjectionWithAttenuation(nowRecEmitIm, attenIm,
							emitRatio, attenRatio);
  Projection ratio(geo);
  for (int i = 0; i < geo.projectionCount; i++) {
    for (int j = 0; j < geo.detectorPixelNum; j++) {
      Real_t tmp = data.getValue(j,i) / (fproj.getValue(j,i) +0.00001);
      ratio.setValue(j, i, tmp);
    }
  }

  //BP
  Reconstruction numerator =op.backwardProjection(ratio);
  Projection flat(geo);
  flat.clear(1.0);
  Reconstruction denominator =op.backwardProjection(flat);

  Reconstruction ret(geo);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      Real_t val =nowRecEmitIm.getValue(x,y)* numerator.getValue(x,y) / denominator.getValue(x,y);
      ret.setValue(x,y, std::max<Real_t>(val, 0.0));
    }
  }

  return std::move(ret);
}



/****************************************
 * OSEM Loop
 */
Reconstruction OSEMLoop(int subsetNum,
			int iteration,
			const Geometry &geo,
			const Reconstruction &nowRecEmitIm,
			const Reconstruction &attenIm,
			double emitRatio,
			double attenRatio,
			const Projection &data)
{
  int start_proj =iteration %subsetNum;
  EmissionProjectionOperator op(geo);

  //FP
  Projection fproj =op.forwardProjectionWithAttenuation(nowRecEmitIm, attenIm,
							emitRatio, attenRatio);
  Projection ratio(geo);
  for (int i = start_proj; i < geo.projectionCount; i+=subsetNum) {
    for (int j = 0; j < geo.detectorPixelNum; j++) {
      Real_t tmp = data.getValue(j,i) / (fproj.getValue(j,i) +0.00001);
      ratio.setValue(j, i, tmp);
    }
  }

  //BP
  Reconstruction numerator =op.backwardProjection(ratio);
  Projection flat(geo);
  flat.clear(1.0);
  Reconstruction denominator =op.backwardProjection(flat);

  Reconstruction ret(geo);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      Real_t val =nowRecEmitIm.getValue(x,y)* numerator.getValue(x,y) / denominator.getValue(x,y);
      ret.setValue(x,y, std::max<Real_t>(val, 0.0));
    }
  }

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
  int maxIteration =40;
  char withInitDir[MBF] ="undefined";
  double init_gain =270000;

  arg.addOption("-output_dir", outputDir, MBF);
  arg.addOption("-emit_ratio", emitRatio);
  arg.addOption("-atten_ratio", attenRatio);
  arg.addOption("-max_iteration", maxIteration);
  arg.addOption("-with_init_dir", withInitDir, MBF);
  arg.addOption("-init_gain", init_gain);

  int retcd = arg.parse(argc, argv);
  if (retcd || argc != 1) {
    printf("Usage: %s\n", argv[0]);
    return -1;
  }


  // general
  char   geometryPath[MBF]           = "../config/geometry_128.ini";
  double noiseIntensity              = 0.25;
  Geometry geo(geometryPath);

  auto petFiles =DirectryOpener(std::string("../Registration/reg_images/pet_test"), 1);
  auto ctFiles =DirectryOpener(std::string("../Registration/reg_images/ct_test"), 1);

  int studyNum =petFiles.getFileNum();
  Real_t tot_cnt_mean =0;
  Real_t tot_cnt_min =99999999999;
  Real_t tot_cnt_max =0;
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

    EmissionProjectionOperator epo(geo);
    Projection emitProj =epo.forwardProjectionWithAttenuation(emitIm, attenIm,
							      emitRatio, attenRatio);
    int seed =i;
    AddCountNoise(geo, &emitProj, seed);
    std::cout << "emit proj:" << std::endl;
    emitProj.stat();
    {
      Real_t tot =emitProj.total_count();
      tot_cnt_mean +=tot /studyNum;
      tot_cnt_min =std::min(tot_cnt_min, tot);
      tot_cnt_max =std::max(tot_cnt_max, tot);
    }


    /****************************************
     * ML-EM
     */
    int iteration =maxIteration +1;
    char algo_str[MBF];
    if(0){
      Reconstruction nowRecEmitIm(geo);
      if(!strcmp(withInitDir,"undefined")){
	sprintf(algo_str, "MLEM");
	InitializePixel(geo, emitProj, &nowRecEmitIm);
      }else{
	sprintf(algo_str, "MLEM_with_init");
	char initNpyName[MBF];
	sprintf(initNpyName, "%s/%s", withInitDir, filename.c_str());
	std::cout << "init npy name:" << initNpyName << std::endl;
	nowRecEmitIm =Npy2ReconstructionFloat32(initNpyName, geo).scale(init_gain);
	nowRecEmitIm.stat();
      }

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =MLEMLoop(geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 10){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 20){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 50){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }


    /****************************************
     * OS-EM
     */
    if(0){
      char   algo_str[MBF];
      int subsetNum =16;
      Reconstruction nowRecEmitIm(geo);
      if(!strcmp(withInitDir,"undefined")){
	sprintf(algo_str, "OSEM");
	InitializePixel(geo, emitProj, &nowRecEmitIm);
      }else{
	sprintf(algo_str, "OSEM_with_init");
	char initNpyName[MBF];
	sprintf(initNpyName, "%s/%s", withInitDir, filename.c_str());
	std::cout << "init npy name:" << initNpyName << std::endl;
	nowRecEmitIm =Npy2ReconstructionFloat32(initNpyName, geo).scale(init_gain);
	nowRecEmitIm.stat();
      }

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =OSEMLoop(subsetNum, iter, geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 10){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 20){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 50){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }



    /****************************************
     * ML-EM with iterative smoothing
     */
    if(0){
      char   algo_str[MBF]           = "MLEM_smooth";
      Real_t sigma =0.45;
      Real_t filtWidth =2;
      Reconstruction nowRecEmitIm(geo);
      InitializePixel(geo, emitProj, &nowRecEmitIm);

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =MLEMLoop(geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.gaussian(sigma, filtWidth);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 10){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 20){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 50){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }


    /****************************************
     * OS-EM with iterative smoothing
     */
    if(1){
      char   algo_str[MBF]           = "OSEM_smooth";
      int subsetNum =16;
      Real_t sigma =0.45;
      Real_t filtWidth =2;
      Reconstruction nowRecEmitIm(geo);
      InitializePixel(geo, emitProj, &nowRecEmitIm);

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =OSEMLoop(subsetNum, iter, geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.gaussian(sigma, filtWidth);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 50){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }


    /****************************************
     * ML-EM with smoothing
     */
    if(0){
      char   algo_str[MBF]           = "MLEM_smooth";
      Real_t sigma =0.15;
      Real_t filtWidth =2;
      Reconstruction nowRecEmitIm(geo);
      InitializePixel(geo, emitProj, &nowRecEmitIm);

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =MLEMLoop(geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 50){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }


    /****************************************
     * OS-EM with smoothing
     */
    if(0){
      char   algo_str[MBF]           = "OSEM_smooth";
      int subsetNum =16;
      Real_t sigma =0.15;
      Real_t filtWidth =2;
      Reconstruction nowRecEmitIm(geo);
      InitializePixel(geo, emitProj, &nowRecEmitIm);

      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRecEmitIm =OSEMLoop(subsetNum, iter, geo, nowRecEmitIm, attenIm, emitRatio, attenRatio, emitProj);
	nowRecEmitIm.stat();
	fprintf(stderr, "\r%s ==%3d-iteration: %f sec.==\n",
		filename.c_str(), iter, (double)(clock() - start) / CLOCKS_PER_SEC );

	char outPath[MBF];
	sprintf(outPath, "./%s/%s/iter%d/%s", outputDir, algo_str, iter, filename.c_str());
	if(iter == 50){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 100){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
	if(iter == 200){ nowRecEmitIm.gaussian(sigma, filtWidth); Reconstruction2Npy(outPath, nowRecEmitIm); }
      }
    }




    /****************************************
     * FBP
     */
    if(0){
      char outPath[MBF];
      sprintf(outPath, "./%s/%s", outputDir, filename.c_str());
      Reconstruction fbpRec =FBP(geo, emitProj);
      std::cout << "FBP rec stat:" << std::endl;
      fbpRec.stat();
      Reconstruction2Npy(outPath, fbpRec);
    }
  }


  std::string output_filename = std::string(outputDir) +std::string("_count.txt");
  std::ofstream outputfile(output_filename.c_str());
  outputfile << "total count mean" << tot_cnt_mean << std::endl;
  outputfile << "total_count_min" << tot_cnt_min << std::endl;
  outputfile << "total_count_max" << tot_cnt_max << std::endl;
  outputfile << std::endl;
  outputfile.close();

}
