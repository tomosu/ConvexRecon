#include "GradientFilter.hpp"
#include "ChambolleProjectionOperator.hpp"
#include <cmath>


Reconstruction ChambolleProjectionOperator(const Real_t lambda, const Real_t step, const Real_t tolerance, const int maxIteration,
					   const Geometry &geo, const Reconstruction &input){

  Reconstruction initial =input;
  Reconstruction lambda_f =input.scale(lambda);

  ReconstructionGradient now_p(geo);
  ReconstructionGradient next_p(geo);
  Real_t maxdiff;

  std::cout << "iteration start" <<std::endl;
  int iter=0;
  do{
    maxdiff =0.0;
    Reconstruction div_p = DivergenceFilter(geo, now_p);
    ReconstructionGradient factor =GradientFilter( geo, div_p.subtract(lambda_f) );
    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	Real2_t term =factor.getValue(x,y).scale(step);
	Real2_t nume =now_p.getValue(x,y) +term;
	Real_t denomi =1.0 +term.norm_L2();
	next_p.setValue(x, y, nume.scale(1.0/denomi)) ;
      }
    }

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	Real2_t xy =next_p.getValue(x,y)-now_p.getValue(x,y);
	maxdiff =std::max(maxdiff, fabs(xy.norm_L2()));
      }
    }
    now_p = next_p;
    iter++;
    std::cout << "iter:" << iter << " diff:" << maxdiff << std::endl;
  }while(maxdiff > tolerance && iter < maxIteration);

  Reconstruction div_p_final =DivergenceFilter(geo, next_p).scale(1.0/lambda);
  Reconstruction ret =initial.subtract(div_p_final);
  return std::move(ret);
}
