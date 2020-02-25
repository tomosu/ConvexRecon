#include "NpyReconConverter.hpp"


Reconstruction Npy2Reconstruction(const char* filename, const Geometry &geo)
{
  cnpy::NpyArray npy_in =cnpy::npy_load(std::string(filename));
  int ny =npy_in.shape[0];
  int nx =npy_in.shape[1];
  double* data =reinterpret_cast<double*>(npy_in.data);
  Reconstruction ret(geo);

  for(int y=0; y<ny; y++){
    for(int x=0; x<nx; x++){
      int idx =x +nx*y;
      ret.setValue(x,y,data[idx]);
    }
  }
  return std::move(ret);
}


Reconstruction Npy2ReconstructionFloat32(const char* filename, const Geometry &geo)
{
  cnpy::NpyArray npy_in =cnpy::npy_load(std::string(filename));
  int ny =npy_in.shape[0];
  int nx =npy_in.shape[1];
  float* data =reinterpret_cast<float*>(npy_in.data);
  Reconstruction ret(geo);

  for(int y=0; y<ny; y++){
    for(int x=0; x<nx; x++){
      int idx =x +nx*y;
      ret.setValue(x,y,double(data[idx]));
    }
  }
  return std::move(ret);
}



void Reconstruction2Npy(const char* filename, Reconstruction in)
{
  unsigned width =in.sizeX;
  unsigned height =in.sizeY;
  std::vector<double> data;

  for(int y=0; y<height; y++){
    for(int x=0; x<width; x++){
      int val =in.getValue(x,y);
      data.push_back((double)val);
    }
  }

  unsigned int shape[] ={height,width};
  cnpy::npy_save(std::string(filename), &data[0], shape, 2, "w");
}
