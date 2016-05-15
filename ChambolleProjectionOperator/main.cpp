/************************************************************
 *
 * Basic 2D-MLEM solver
 *
 ************************************************************/
#include "Geometry.hpp"
#include "FP.hpp"
#include "ProjectionOperator.hpp"
#include "ChambolleProjectionOperator.hpp"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <random>

#define MBF 256


/****************************************
 * initialize pixel
 */
inline void InitializeImage(const Geometry &geo,
			    Reconstruction *recon)
{

  //color bar
  recon->clear(1.0);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      Real_t val = (x < geo.volumeSize.x/3) ? 0.6 : 0.3;
      val = (x > 2*geo.volumeSize.x/3) ? 1.0 : val;
      recon->setValue(x,y,val);
    }
  }

  //add noise
  std::mt19937 mt(0);
  std::uniform_real_distribution<double> genrand1(0.0, 1.0);
  std::uniform_int_distribution<int> genrand_x(0, geo.volumeSize.x-1);
  std::uniform_int_distribution<int> genrand_y(0, geo.volumeSize.y-1);

  for(int i=0; i<40000; i++){
    Real_t noise   = sqrtf(-2 * logf(genrand1(mt)) ) * cosf(2 * M_PI * genrand1(mt)); // BoxMuller
    noise *= 0.2;

    int x =genrand_x(mt);
    int y =genrand_y(mt);
    Real_t val =recon->getValue(x,y) +noise;
    recon->setValue(x,y,val);
  }
}



/****************************************
 * output recon
 */
void OutputBinaryRecon(int   iteration,
		       const Geometry &geo,
		       const Reconstruction& rec,
                       char  *outputDir,
                       int   convertRealToUnsignedShort)
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
    sprintf(filename, "%s/output_%02d_s16.bin", outputDir, iteration);
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
    sprintf(filename, "%s/output_%02d_f32.bin", outputDir, iteration);
    if (FILE * fp = fopen(filename, "wb+")) {
      fwrite(pixel, sizeof(float), geo.volumeSize.x * geo.volumeSize.y, fp);
      fclose(fp);
    }
    delete[] pixel;

  }
}


/****************************************
 * main
 */
int main(int argc, char **argv)
{
  char   phantomPath[MBF]            = "../config/phantom.ini";
  char   geometryPath[MBF]           = "../config/geometry.ini";
  char   outputDir[MBF]              = "./output";

  Geometry geo(geometryPath);
  Reconstruction initial(geo);

  InitializeImage(geo, &initial);
  std::cout << "init done" << std::endl;

  OutputBinaryRecon(0, geo, initial, outputDir, 0);
  std::cout << "output initial" << std::endl;

  Real_t lambda=10.0;
  Real_t step=0.25;
  Real_t tolerance=0.001;
  int maxIteration=1000;

  Reconstruction result =ChambolleProjectionOperator(lambda, step, tolerance, maxIteration, geo, initial);
  OutputBinaryRecon(1, geo, result, outputDir, 0);
}
