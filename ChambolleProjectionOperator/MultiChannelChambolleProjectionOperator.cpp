#include "GradientFilter.hpp"
#include "ChambolleProjectionOperator.hpp"
#include <cmath>


std::vector<Reconstruction> MultiChannelChambolleProjectionOperator(const Real_t lambda, const Real_t step, const Real_t tolerance, const int maxIteration,
								    const Geometry &geo, const std::vector<Reconstruction> &input){

  std::vector<Reconstruction> initial =input;
  std::vector<Reconstruction> lambda_f =input;
  for(auto &l : lambda_f){
    l =l.scale(lambda);
  }

  std::vector<ReconstructionGradient> now_p, next_p;
  for(int i=0; i<input.size(); i++){
    now_p.push_back(ReconstructionGradient(geo));
    next_p.push_back(ReconstructionGradient(geo));
  }
  Real_t maxdiff;

  std::cout << "iteration start" <<std::endl;
  int iter=0;
  do{
    maxdiff =0.0;

    std::vector<ReconstructionGradient> factor;
    std::vector<Reconstruction> div_p;
    for(int i=0; i<input.size(); i++){
      div_p.push_back( DivergenceFilter( geo, now_p[i]) );
      factor.push_back( GradientFilter( geo, div_p[i].subtract(lambda_f[i])) );
    }

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	
	//calc all channel norm
	Real_t termNorm=0.0;
	for(int i=0; i<input.size(); i++){
	  Real2_t term =factor[i].getValue(x,y);
	  termNorm +=term.norm_L2();
	}
	termNorm =sqrt(termNorm) *step;

	//update equation
	Real_t denomi =termNorm +1.0;
	for(int i=0; i<input.size(); i++){
	  Real2_t term =factor[i].getValue(x,y).scale(step);
	  Real2_t nume =now_p[i].getValue(x,y) +term;
	  next_p[i].setValue(x, y, nume.scale(1.0/denomi)) ;
	}
      }
    }

    for(int y=0; y<geo.volumeSize.y; y++){
      for(int x=0; x<geo.volumeSize.x; x++){
	for(int i=0; i<input.size(); i++){
	  Real2_t xy =next_p[i].getValue(x,y)-now_p[i].getValue(x,y);
	  maxdiff =std::max(maxdiff, fabs(xy.norm_L2()));
	}
      }
    }

    now_p = next_p;
    iter++;
    std::cout << "iter:" << iter << " diff:" << maxdiff << std::endl;
  } while(maxdiff > tolerance && iter < maxIteration);

  std::vector<Reconstruction> div_p_final;
  for(int i=0; i<input.size(); i++){
    div_p_final.push_back( DivergenceFilter(geo, next_p[i]).scale(1.0/lambda) );
  }

  std::vector<Reconstruction> ret;
  for(int i=0; i<input.size(); i++){
    ret.push_back( initial[i].subtract(div_p_final[i]) );
  }

  return std::move(ret);
}
