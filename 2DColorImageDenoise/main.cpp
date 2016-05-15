/************************************************************
 *
 * Basic 2D-MLEM solver
 *
 ************************************************************/
#include "Geometry.hpp"
#include "PngReconConverter.hpp"
#include "FP.hpp"
#include "ProjectionOperator.hpp"
#include "MultiChannelChambolleProjectionOperator.hpp"
#include "lodepng.h"

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
inline void AddNoise(const Geometry &geo,
		     Reconstruction *recon,
		     int seed)
{
  //add noise
  std::mt19937 mt(seed);
  std::uniform_real_distribution<double> genrand1(0.0, 1.0);
  std::uniform_int_distribution<int> genrand_x(0, geo.volumeSize.x-1);
  std::uniform_int_distribution<int> genrand_y(0, geo.volumeSize.y-1);

  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      //  for(int i=0; i<800000; i++){
      Real_t noise   = sqrtf(-2 * logf(genrand1(mt)) ) * cosf(2 * M_PI * genrand1(mt)); // BoxMuller
      noise *= 0.1;
      
      //int x =genrand_x(mt);
      //int y =genrand_y(mt);
      
      Real_t val =recon->getValue(x,y) +noise;
      recon->setValue(x,y,val);
    }
  }

}



/****************************************
 * main
 */
int main(int argc, char **argv)
{
  char   geometryPath[MBF]           = "../config/geometry_512.ini";
  char   inputPath[MBF]              = "./in.png";
  char   initialPath[MBF]              = "./initial.png";
  char   outputPath[MBF]              = "./out.png";

  Geometry geo(geometryPath);
  std::vector<Reconstruction> initial =Png2Reconstruction(inputPath, geo);

  for(int i=0; i<4; i++){
    AddNoise(geo, &initial[i], i);
  }
  Reconstruction2Png(initialPath, geo, initial);
  std::cout << "init done" << std::endl;

  //Real_t lambda=10.0;
  //Real_t lambda=3.0;
  Real_t lambda=1.0;
  //Real_t lambda=0.5;
  Real_t step=0.25;
  Real_t tolerance=0.0001;
  int maxIteration=100;

  std::vector<Reconstruction> result =MultiChannelChambolleProjectionOperator(lambda, step, tolerance, maxIteration, geo, initial);
  Reconstruction2Png(outputPath, geo, result);
}
