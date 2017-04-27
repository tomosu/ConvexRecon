#pragma once

#include "ProjectionOperator.hpp"
#include <cstdint>
#include <string>

class DicomSliceReader{

public:
  DicomSliceReader(const char *filepath, int sizeX_in, int sizeY_in)
    : sizeX(sizeX_in), sizeY(sizeY_in)
  {
    int size =sizeX *sizeY;
    rawPixelMax =INT16_MIN;
    rawPixelMin =INT16_MAX;
    rawPixels.resize(size);

    FILE* fp =fopen(filepath, "rb");
    if(fp != NULL){
      fseek(fp, -size*2, SEEK_END);
      fread(&rawPixels[0], sizeof(int16_t), size, fp);
      fclose(fp);
    }else{
      std::cout << "invalid file " << std::string(filepath) << std::endl;
      exit(0);
    }

    normalizedPixels.resize(size);
    Normalize();
  }

  ~DicomSliceReader(){  }

  Reconstruction generateReconstruction(const Geometry &geoIn){
    assert(sizeX == geoIn.volumeSize.x);
    assert(sizeY == geoIn.volumeSize.y);
    Reconstruction ret(geoIn);
    for(int y=0; y<sizeY; y++){
      for(int x=0; x<sizeX; x++){
	int idx =x+sizeX*y;
	ret.setValue(x, y, normalizedPixels[idx]);
      }
    }
    return std::move(ret);
  }


private:
  void Normalize(){
    for(auto p : rawPixels){
      rawPixelMax =std::max(p, rawPixelMax);
      rawPixelMin =std::min(p, rawPixelMin);
    }
    for(int i=0; i<rawPixels.size(); i++){
      normalizedPixels[i] =(Real_t)(rawPixels[i] -rawPixelMin);
    }
  }

private:
  int sizeX;
  int sizeY;
  std::vector<int16_t> rawPixels;
  std::vector<Real_t> normalizedPixels;
  int16_t rawPixelMin;
  int16_t rawPixelMax;
};
