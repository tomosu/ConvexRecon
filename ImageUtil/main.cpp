#include "ProjectionOperator.hpp"
#include "PngReconConverter.hpp"
#include "NpyReconConverter.hpp"
#include "lodepng.h"

#define MBF 256

#if 0
int main(){
  char   geometryPath[MBF]        = "../config/geometry_512.ini";
  Geometry geo(geometryPath);

  char   inputPath[MBF]           = "./in.png";
  char   outputPath[MBF]          = "./out.png";

  std::vector<Reconstruction> image =Png2Reconstruction(inputPath, geo);
  std::cout << "load" << std::endl;
  Reconstruction2Png(outputPath, geo, image);
  std::cout << "save" << std::endl;
  return 0;
}
#endif


#if 1
int main(){
  char   geometryPath[MBF]        = "../config/geometry_128.ini";
  Geometry geo(geometryPath);

  char   inputPath[MBF]           = "./in.npy";
  char   outputPath[MBF]          = "./out.npy";

  Reconstruction rec =Npy2Reconstruction(inputPath, geo);
  std::cout << "load" << std::endl;

  Reconstruction2Npy(outputPath, rec);
  std::cout << "save" << std::endl;

  Reconstruction rec2 =Npy2Reconstruction(outputPath, geo);
  std::cout << "load" << std::endl;
  return 0;
}
#endif
