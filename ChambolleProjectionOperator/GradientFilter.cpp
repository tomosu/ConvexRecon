#include "Geometry.hpp"
#include "GradientFilter.hpp"
#include "ProjectionOperator.hpp"
#include <cmath>


//adjoint operator of divergence filter
ReconstructionGradient GradientFilter(const Geometry &geo,
				      const Reconstruction &in){
  ReconstructionGradient ret(geo);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){

      Real_t diff_x = (x == geo.volumeSize.x-1) ? 0.0 : in.getValue(x+1,y) -in.getValue(x,y);
      ret.setDxValue(x, y, diff_x);

      Real_t diff_y = (y == geo.volumeSize.y-1) ? 0.0 : in.getValue(x,y+1) -in.getValue(x,y);
      ret.setDyValue(x, y, diff_y);
    }
  }
  return std::move(ret);
}


//adjoint operator of gradient filter
Reconstruction DivergenceFilter(const Geometry &geo,
				const ReconstructionGradient &in){

  Reconstruction ret(geo);
  for(int y=0; y<geo.volumeSize.y; y++){
    for(int x=0; x<geo.volumeSize.x; x++){
      
      Real_t diff_x = (x == 0) ? in.getDxValue(x,y) : in.getDxValue(x,y) -in.getDxValue(x-1,y);
      diff_x = (x == geo.volumeSize.x-1 && x > 0 ) ? -in.getDxValue(x-1,y) : diff_x;
      
      Real_t diff_y = (y == 0) ? in.getDyValue(x,y) : in.getDyValue(x,y) -in.getDyValue(x,y-1);
      diff_y = (y == geo.volumeSize.y-1 && y > 0 )? -in.getDyValue(x,y-1) : diff_y;
      
      ret.setValue(x,y, diff_x +diff_y);
    }
  }

  return std::move(ret);
}
