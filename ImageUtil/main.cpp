#include "ProjectionOperator.hpp"
#include "PngReconConverter.hpp"
#include "lodepng.h"

#define MBF 256

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
