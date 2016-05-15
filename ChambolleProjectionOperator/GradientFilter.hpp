#pragma once

#include "Geometry.hpp"
#include "ProjectionOperator.hpp"
#include <cmath>


struct ReconstructionGradient{

  Reconstruction dx;
  Reconstruction dy;
  
  ReconstructionGradient(){}

  ReconstructionGradient(const Geometry &geo)
    : dx{Reconstruction(geo)}, dy{Reconstruction(geo)} {}


  ReconstructionGradient(const ReconstructionGradient& p){
    dx =p.dx;
    dy =p.dy;
  }

  ReconstructionGradient(ReconstructionGradient&& p){
    dx =std::move(p.dx);
    dy =std::move(p.dy);
  }

  ReconstructionGradient& operator=(const ReconstructionGradient& p){
    dx =p.dx;
    dy =p.dy;
    return *this;
  }

  ReconstructionGradient& operator=(ReconstructionGradient&& p){
    dx =std::move(p.dx);
    dy =std::move(p.dy);
    return *this;
  }

  void setValue(const int x, const int y, const Real2_t &val){dx.setValue(x,y,val.x); dy.setValue(x,y,val.y);}
  void setDxValue(const int x, const int y, const Real_t val){dx.setValue(x,y,val);}
  void setDyValue(const int x, const int y, const Real_t val){dy.setValue(x,y,val);}

  Real2_t getValue(const int x, const int y) const { return Real2_t( dx.getValue(x,y), dy.getValue(x,y) );}
  Real_t getDxValue(const int x, const int y) const { return dx.getValue(x,y);}
  Real_t getDyValue(const int x, const int y) const { return dy.getValue(x,y);}
};

//adjoint operator of divergence filter
ReconstructionGradient GradientFilter(const Geometry &geo,
				      const Reconstruction &in);

//adjoint operator of gradient filter
Reconstruction DivergenceFilter(const Geometry &geo,
				const ReconstructionGradient &in);
