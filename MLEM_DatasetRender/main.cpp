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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <algorithm>
#include <complex>

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
			const Reconstruction &nowRec,
			const Projection &data)
{
  ProjectionOperator op(geo);

  //FP
  Projection fproj =op.forwardProjection(nowRec);
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
      Real_t val =nowRec.getValue(x,y)* numerator.getValue(x,y) / denominator.getValue(x,y);
      ret.setValue(x,y, std::max<Real_t>(val, 0.0));
    }
  }

  return std::move(ret);
}


/****************************************
 * output recon
 */
void OutputBinaryRecon(const Geometry &geo,
		       const Reconstruction& rec,
                       char   *outputDir,
		       char   *outputFileName,
                       bool convertRealToUnsignedShort)
{
  if (convertRealToUnsignedShort) {

    unsigned short *pixel = new unsigned short[geo.volumeSize.x * geo.volumeSize.y]();

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	int i =x +geo.volumeSize.x*y;
	pixel[i] = (unsigned short)floor(std::max<Real_t>(rec.getValue(x,y), 0.0) * 100000);
      }
    }

    char filename[256];
    sprintf(filename, "%s/%s.bin", outputDir, outputFileName);
    if (FILE * fp = fopen(filename, "wb+")) {
      fwrite(pixel, sizeof(unsigned short), geo.volumeSize.x * geo.volumeSize.y, fp);
      fclose(fp);
    }
    delete[] pixel;

  } else {

    Real_t *pixel = new Real_t[geo.volumeSize.x * geo.volumeSize.y]();

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	int i =x +geo.volumeSize.x*y;
	pixel[i] = (Real_t)(std::max<Real_t>(rec.getValue(x,y), 0.0));
      }
    }

    char filename[256];
    sprintf(filename, "%s/%s.bin", outputDir, outputFileName);
    if (FILE * fp = fopen(filename, "wb+")) {
      fwrite(pixel, sizeof(Real_t), geo.volumeSize.x * geo.volumeSize.y, fp);
      fclose(fp);
    }
    delete[] pixel;

  }
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

  std::cout<<"OpenMP : Enabled (Max # of threads = "<<omp_get_max_threads()<<")"<< std::endl;
  // general
  char   geometryPath[MBF]           = "../config/geometry_512.ini";
  double noiseIntensity              = 0.25;
  double digitizedAir               = 16383;
  Geometry geo(geometryPath);

  // filepaths
  int dcmLevel =4;
  DicomDirectryOpener dopener(std::string("../DICOM/pancreas/DOI"), dcmLevel);
  int dataNum =dopener.getFileNum();
  std::cout << "file num:" << dataNum << std::endl;

  for(int i=0; i<dataNum; i++){
    std::cout << "current image index:" << i << std::endl;

    /****************************************
     * set answer
     */
    std::string dicomPath = dopener.getFilePath(i);
    DicomSliceReader dcmReader(dicomPath.c_str(),
			       geo.volumeSize.x,
			       geo.volumeSize.y);
    Reconstruction answer =dcmReader.generateReconstruction(geo);
    answer =answer.scale(0.00001); //tekitou
    answer.stat();
    {
      char answerOutputDir[MBF];
      sprintf(answerOutputDir, "./results/answer");
      char filename[MBF];
      sprintf(filename, "%s", CreateOmittedFilename(dicomPath).c_str());
      OutputBinaryRecon(geo, answer, answerOutputDir, filename, false);
    }


    //clean or noise
    for(int j=0; j<2; j++){
      /****************************************
       * generate projection
       */
      bool isNoisy =(j%2 == 1);
      char outputDir[MBF];
      sprintf(outputDir, "%s", isNoisy ? "./results/noise_output" : "./results/output");

      ProjectionOperator dataop(geo);
      Projection data =dataop.forwardProjection(answer);

      if(isNoisy){
	unsigned int seed =i;
	AddNoise(geo, noiseIntensity, digitizedAir, &data, seed);
      }


      /****************************************
       * FBP
       */
      {
	char fbpOutputDir[MBF];
	sprintf(fbpOutputDir, "%s", isNoisy ? "./results/noise_output_fbp" : "./results/output_fbp");
	char filename[MBF];
	sprintf(filename, "%s", CreateOmittedFilename(dicomPath).c_str());
	Reconstruction fbpRec =FBP(geo, data);
	std::cout << "FBP rec stat:" << std::endl;
	fbpRec.stat();
	OutputBinaryRecon(geo, fbpRec, fbpOutputDir, filename, false);
      }


      /****************************************
       * fill pixels with appropriate value
       */
      Reconstruction nowRec(geo);
      InitializePixel(geo, data, &nowRec);


      /****************************************
       * MLEM reconstruction
       */
      int iteration =201;
      char filename[MBF];
      sprintf(filename, "%s", CreateOmittedFilename(dicomPath).c_str());
      for (int iter = 0; iter < iteration; iter++) {
	clock_t start = clock();
	nowRec = MLEMLoop(geo, nowRec, data);
	nowRec.stat();
	char cleanOrNoisy[MBF];
	sprintf(cleanOrNoisy, "%s", isNoisy ? "noisy" : "clean");
	fprintf(stderr, "\r%s %s ==%3d-iteration: %f sec.==\n",
		filename, cleanOrNoisy, iter,
		(double)(clock() - start) / CLOCKS_PER_SEC );


	/// output in specific iteration
	char iterDir[MBF];
	sprintf(iterDir, "%s/iter%d", outputDir, iter);
	if(iter == 10){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}

	if(iter == 25){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}

	if(iter == 50){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}

	if(iter == 100){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}

	if(iter == 150){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}

	if(iter == 200){
	  OutputBinaryRecon(geo, nowRec, iterDir, filename, false);
	}
      }
    }
  }
}
