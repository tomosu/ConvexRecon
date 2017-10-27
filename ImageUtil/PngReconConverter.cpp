#include "PngReconConverter.hpp"

std::vector<Reconstruction> Png2Reconstruction(const char* filename, const Geometry &geo)
{
  std::vector<unsigned char> image; //the raw pixels
  unsigned width, height;

  //decode
  unsigned error = lodepng::decode(image, width, height, filename);

  //if there's an error, display it
  if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
  if( !(geo.volumeSize.x == width && geo.volumeSize.y == height) ){
    std::cout << "invalid size" << std::endl;
    exit(0);
  }

  std::vector<Reconstruction> ret;
  for(int i=0; i<4; i++){
    ret.push_back(Reconstruction(geo));
  }

  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      for(int i=0; i<4; i++){
	Real_t val =(Real_t)image[4*(x +width*y) +i]/(Real_t)255;
	ret[i].setValue(x,y,val);
      }
    }
  }

  return std::move(ret);
}



void Reconstruction2Png(const char* filename, const Geometry &geo, std::vector<Reconstruction> in)
{
  std::vector<unsigned char> image; //the raw pixels

  unsigned width =geo.volumeSize.x;
  unsigned height =geo.volumeSize.y;

  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      for(int i=0; i<4; i++){
	int ival =(in[i].getValue(x,y)*255);
	unsigned char uval =(unsigned char)std::max( std::min(255, ival), 0);
	image.push_back(uval);
      }
    }
  }

  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

}
