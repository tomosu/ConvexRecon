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
      ret.setValue(x,y, std::max(val, 0.0));
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
                       int     convertRealToUnsignedShort)
{
  if (convertRealToUnsignedShort) {

    unsigned short *pixel = new unsigned short[geo.volumeSize.x * geo.volumeSize.y]();

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	int i =x +geo.volumeSize.x*y;
	pixel[i] = (unsigned short)floor(std::max(rec.getValue(x,y), 0.0) * 100000);
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

    float *pixel = new float[geo.volumeSize.x * geo.volumeSize.y]();

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	int i =x +geo.volumeSize.x*y;
	pixel[i] = (float)(std::max(rec.getValue(x,y), 0.0));
      }
    }

    char filename[256];
    sprintf(filename, "%s/%s.bin", outputDir, outputFileName);
    if (FILE * fp = fopen(filename, "wb+")) {
      fwrite(pixel, sizeof(float), geo.volumeSize.x * geo.volumeSize.y, fp);
      fclose(fp);
    }
    delete[] pixel;

  }
}


/****************************************
 *
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
  double noiseIntensity              = 0.2;
  double digitizedAir               = 16383;
  Geometry geo(geometryPath);

  // filepaths
  int dcmLevel =4;
  DicomDirectryOpener dopener(std::string("../DICOM/pancreas/DOI"), dcmLevel);
  int dataNum =dopener.getFileNum();
  std::cout << "file num:" << dataNum << std::endl;
  for(int i=130; i<dataNum; i++){
    std::cout << "current image index:" << i << std::endl;
    //clean or noise
    for(int j=0; j<2; j++){
      /****************************************
       * generate projection
       */
      bool isNoisy =(j%2 == 1);
      char outputDir[MBF];
      sprintf(outputDir, "%s", isNoisy ? "./noise_output" : "./output");
      std::string dicomPath = dopener.getFilePath(i);
      DicomSliceReader dcmReader(dicomPath.c_str(),
				 geo.volumeSize.x,
				 geo.volumeSize.y);
      Reconstruction answer =dcmReader.generateReconstruction(geo);
      answer =answer.scale(0.00001); //tekitou
      answer.stat();

      ProjectionOperator dataop(geo);
      Projection data =dataop.forwardProjection(answer);
      if(isNoisy){
	AddNoise(geo, noiseIntensity, digitizedAir, &data);
      }

      /****************************************
       * fill pixels with appropriate value
       */
      Reconstruction nowRec(geo);
      InitializePixel(geo, data, &nowRec);

      /****************************************
       * iterate reconstruction
       */
      int iteration =200;
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
      }
      OutputBinaryRecon(geo, nowRec, outputDir, filename, 0);
    }
  }
}
