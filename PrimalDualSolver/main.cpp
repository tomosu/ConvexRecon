/************************************************************
 *
 * Primal Dual Method
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
#include "GradientFilter.hpp"
#include "ChambolleProjectionOperator.hpp"

#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
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
 * Recon Loop
 */
std::tuple<Reconstruction,Projection> APPLoop(const Geometry &geo,
					      const Real_t primalStep,
					      const Real_t dualStep,
					      const Projection &nowL,
					      const Reconstruction &nowRec,
					      const Projection &data)
{
  ProjectionOperator op(geo);

  //step 2-1
  Projection fproj =op.forwardProjection(nowRec);
  Projection residual =fproj.subtract(data);
  Projection tempL =nowL.add( residual.scale(dualStep) );

  //step 2-2
  Reconstruction bp =op.backwardProjection(tempL);
  Reconstruction proxIn =nowRec.subtract( bp.scale(primalStep) );
  Real_t lambda=2000.0;
  Real_t cStep=0.25;
  Real_t tol=0.00001;
  Real_t maxIter=100;
  Reconstruction nextRec =ChambolleProjectionOperator(lambda, cStep, tol, maxIter, geo, proxIn);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      Real_t val =std::max(0.0 ,nextRec.getValue(x,y));
      nextRec.setValue(x,y,val);
    }
  }

  //step 2-3
  Projection nextFproj =op.forwardProjection(nextRec);
  Projection nextResidual =nextFproj.subtract(data);
  Projection nextL =nowL.add( nextResidual.scale(dualStep) );
  std::tuple<Reconstruction,Projection> ret =std::make_tuple(nextRec, nextL);
  return std::move(ret);
}


/****************************************
 * output recon
 */
void OutputBinaryRecon(int     iteration,
		       const Geometry &geo,
		       const Reconstruction& rec,
                       char   *outputDir,
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

// void OutputTextRecon(Image_t rec,
//                      Real_t  *pixelValue,
//                      char   *outputDir,
//                      char   *outputFile)
// {
//   char fn[MBF];
//   sprintf(fn, "%s/%s", outputDir, outputFile);

//   FILE *fp = fopen(fn, "wb+");
//   if (fp == NULL) {
//     printf("recon file error\n");
//     exit(0);
//   }
//   for (int y = 0; y < geo.volumeSize.y; y++) {
//     for (int x = 0; x < geo.volumeSize.x; x++) {
//       Real2_t pixelCenter;
//       pixelCenter.x = rec.center.x + (x - geo.volumeSize.x / 2) * rec.pitch.x + rec.pitch.x / 2.0;
//       pixelCenter.y = rec.center.y + (y - geo.volumeSize.y / 2) * rec.pitch.y + rec.pitch.y / 2.0;

//       int pixelIndex = y * geo.volumeSize.x + x;
//       fprintf(fp, "%f %f %f\n", pixelCenter.x, pixelCenter.y, pixelValue[pixelIndex]);
//     }
//     fprintf(fp, "\n");
//   }
//   fclose(fp);
// }

/***************************************
 * output phantom
 */

// #define INSIDE  1
// #define OUTSIDE 0
// int IsInside(Phantom_t phantom,
//              Real2_t  pixelPos)
// {
//   Real_t rad2  = phantom.rad * phantom.rad;
//   Real_t dist2 = (phantom.center.x - pixelPos.x) * (phantom.center.x - pixelPos.x) +
//     (phantom.center.y - pixelPos.y) * (phantom.center.y - pixelPos.y);

//   if (dist2 < rad2) {
//     return INSIDE;
//   } else {
//     return OUTSIDE;
//   }
// }

// void OutputTextPhantom(Image_t rec,
//                        char   *phantomPath,
//                        char   *outputDir,
//                        char   *outputFile)
// {
//   char fn[MBF];
//   sprintf(fn, "%s/%s", outputDir, outputFile);

//   FILE *fp = fopen(fn, "wb+");
//   if (fp == NULL) {
//     printf("phantom file error\n");
//     exit(0);
//   }

//   int        phantomNum = count_body(phantomPath);
//   Phantom_t *phantoms   = setupPhantom(phantomPath);

//   for (int y = 0; y < geo.volumeSize.y; y++) {
//     for (int x = 0; x < geo.volumeSize.x; x++) {
//       Real2_t pixelPosition;
//       pixelPosition.x = rec.center.x + (x - geo.volumeSize.x / 2) * rec.pitch.x + rec.pitch.x / 4.0;
//       pixelPosition.y = rec.center.y + (y - geo.volumeSize.y / 2) * rec.pitch.y + rec.pitch.y / 4.0;

//       Real2_t pixelCellCenter[4];
//       for (int i = 0; i < 4; i++) {
// 	pixelCellCenter[i].x = pixelPosition.x + ((i       % 2) == 1) * rec.pitch.x / 2.0;
// 	pixelCellCenter[i].y = pixelPosition.y + (((i / 2) % 2) == 1) * rec.pitch.y / 2.0;
//       }

//       Real_t val = 0.0;
//       for (int i = 0; i < phantomNum; i++) {
// 	// normal signal

// 	if (phantoms[i].isBase == NOT_BASE) {
// 	  for (int j = 0; j < 4; j++) {
// 	    if (IsInside(phantoms[i],
// 			 pixelCellCenter[j]) == INSIDE) {
// 	      val += phantoms[i].mu * 0.25;
// 	    }
// 	  }
// 	}
//       }

//       for (int i = 0; i < phantomNum; i++) {
// 	// base signal

// 	if (val == 0.0) {
// 	  if (phantoms[i].isBase == BASE) {
// 	    for (int j = 0; j < 4; j++) {
// 	      if (IsInside(phantoms[i],
// 			   pixelCellCenter[j]) == INSIDE) {
// 		val += phantoms[i].mu * 0.25;
// 	      }
// 	    }
// 	  }
// 	}
//       }

//       Real2_t pixelCenter;
//       pixelCenter.x = rec.center.x + (x - geo.volumeSize.x / 2) * rec.pitch.x + rec.pitch.x / 2.0;
//       pixelCenter.y = rec.center.y + (y - geo.volumeSize.y / 2) * rec.pitch.y + rec.pitch.y / 2.0;

//       fprintf(fp, "%f %f %f\n", pixelCenter.x, pixelCenter.y, val);
//     }
//     fprintf(fp, "\n");
//   }
//   fclose(fp);
// }


/****************************************
 * main
 */
int main(int argc, char **argv)
{
  char   phantomPath[MBF]            = "../config/phantom.ini";
  char   geometryPath[MBF]           = "../config/geometry.ini";
  char   outputDir[MBF]              = "./output";
  double noiseIntensity              = 1.0;
  double digitizedAir               = 16383;

  Geometry geo(geometryPath);
  /****************************************
   * generate projection
   */
  Projection data =CalcProjectionValue(geo, phantomPath);
  AddNoise(geo, noiseIntensity, digitizedAir, &data);
  /****************************************
   * fill pixels with appropriate value
   */
  Reconstruction nowRec(geo);
  InitializePixel(geo, data, &nowRec);
  Projection nowL(geo);

  /****************************************
   * iterate reconstruction
   */
  int iteration =1000;
  for (int i = 0; i < iteration; i++) {
    Real_t primalStep =0.1/sqrt((Real_t)i/2.0+1.0);
    Real_t dualStep =0.1/sqrt((Real_t)i/2.0+1.0);
    clock_t start = clock();
    auto ret =APPLoop(geo, primalStep, dualStep, nowL, nowRec, data);
    nowRec =std::get<0>(ret);
    nowL =std::get<1>(ret);
    nowRec.stat();
    OutputBinaryRecon(i, geo, nowRec, outputDir, 0);
    fprintf(stderr, "\r=========%3d-iteration: %f sec.===========\n", i, (double)(clock() - start) / CLOCKS_PER_SEC );
  }

  /****************************************
   * convert to text
   */
  //OutputTextPhantom(rec, phantomPath, outputDir, phantomOutText);
  //OutputTextRecon(rec, pixelValue, outputDir, reconOutText);
}
