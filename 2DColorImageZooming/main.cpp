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
inline std::vector<Reconstruction> Shrink(const int N,
					  const Geometry &geo,
					  const std::vector<Reconstruction> recon)
{
  std::vector<Reconstruction> ret;
  for(int i=0; i<recon.size(); i++){
    ret.push_back(Reconstruction(geo));
    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	Real_t val=0.0;
	int idx_x =x-x %N;
	int idx_y =y-y %N;
	for(int dy=0; dy<N; dy++){
	  for(int dx=0; dx<N; dx++){
	    val +=recon[i].getValue(idx_x +dx, idx_y +dy);
	  }
	}
	val /=Real_t(N*N);
	ret[i].setValue(x,y,val);
      }
    }
  }
  return std::move(ret);
}



/****************************************
 * main
 */
int main(int argc, char **argv)
{
  char   geometryPath[MBF]           = "../config/geometry_512.ini";
  char   inputPath[MBF]              = "./mono_image_Lena512.png";
  char   initialPath[MBF]              = "./initial.png";
  char   outputPath[MBF]              = "./out.png";

  Geometry geo(geometryPath);
  std::vector<Reconstruction> in =Png2Reconstruction(inputPath, geo);
  std::vector<Reconstruction> initial =Shrink(4, geo, in);
  
  Reconstruction2Png(initialPath, geo, initial);
  std::cout << "init done" << std::endl;

  Real_t lambda=300.0;
  //Real_t lambda=30.0;
  //Real_t lambda=10.0;
  //Real_t lambda=3.0;
  //Real_t lambda=1.0;  
  Real_t step=0.25;
  Real_t tolerance=0.001;
  int maxIteration=1000;

  std::vector<Reconstruction> u_now, w_now, g_w_now;
  w_now.resize(4);
  g_w_now.resize(4);
  u_now.resize(4);
  for(int i=0; i<4; i++){
    w_now[i]=Reconstruction(geo);
    g_w_now[i]=Reconstruction(geo);
    u_now[i]=Reconstruction(geo);
  }

  //iteration
  for(int i=0; i<100; i++){
    for(int j=0; j<4; j++){ g_w_now[j] =w_now[j].add(initial[j]); }
    u_now =MultiChannelChambolleProjectionOperator(lambda, step, tolerance, maxIteration, geo, g_w_now);
    for(int j=0; j<4; j++){ w_now[j] =u_now[j].subtract(initial[j]); }
    sprintf(outputPath, "./lenna_4x_output/out_loop%i.png", i);
    Reconstruction2Png(outputPath, geo, u_now);
  }

}
