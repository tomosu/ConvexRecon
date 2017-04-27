#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "RealTypes.hpp"
#include "FP.hpp"
#include "Geometry.hpp"


struct Projection{

  std::vector<Real_t> value;
  int detectorPixelNum;
  int projectionCount;

  Projection(){}

  Projection(const Geometry &geoIn){
    detectorPixelNum =geoIn.detectorPixelNum;
    projectionCount =geoIn.projectionCount;
    value.resize(geoIn.detectorPixelNum *geoIn.projectionCount);
    for(auto &v : value ){
      v=0.0;
    }

  }


  Projection(const Projection& p){
    detectorPixelNum =p.detectorPixelNum;
    projectionCount =p.projectionCount;
    value.resize(p.value.size());
    for(int i=0; i<p.value.size(); i++){
      value[i]=p.value[i];
    }
  }

  Projection(Projection&& p){
    detectorPixelNum =p.detectorPixelNum;
    projectionCount =p.projectionCount;
    value =std::move(p.value);
  }

  const Projection& operator=(const Projection& p) {
    detectorPixelNum =p.detectorPixelNum;
    projectionCount =p.projectionCount;
    value.resize(p.value.size());
    for(int i=0; i<p.value.size(); i++){
      value[i]=p.value[i];
    }
    return *this;
  }

  Projection& operator=(Projection&& p){
    detectorPixelNum =p.detectorPixelNum;
    projectionCount =p.projectionCount;
    value =std::move(p.value);
    return *this;
  }

  Projection scale(const Real_t scalar) const {
    Projection ret=*this;
    for(auto &v : ret.value ){
      v *=scalar;
    }
    return std::move(ret);
  }

  Projection subtract(const Projection &r) const {
    Projection ret=*this;
    assert(ret.value.size()==r.value.size());
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]-=r.value[i];
    }
    return std::move(ret);
  }

  Projection add(const Projection &r) const {
    Projection ret=*this;
    assert(ret.value.size()==r.value.size());
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]+=r.value[i];
    }
    return std::move(ret);
  }


  Real_t operator[](int i){
    return value[i];
  }

  Real_t getValue(int detectorPos, int projection) const {
    return value[detectorPos +projection *detectorPixelNum];
  }

  void setValue(const int detectorPos, const int projection, const Real_t val){
    value[detectorPos +projection *detectorPixelNum] =val;
  }

  void clear(double val=0.0){
    for(auto &v : value ){
      v=val;
    }
  }

};



struct Reconstruction{

  std::vector<Real_t> value;
  int sizeX, sizeY;

  Reconstruction(const Geometry &geoIn){
    sizeX =geoIn.volumeSize.x;
    sizeY =geoIn.volumeSize.y;
    value.resize(geoIn.volumeSize.x * geoIn.volumeSize.y);
    for(auto &v : value ){
      v=0.0;
    }
  }

  Reconstruction(){}

  Reconstruction(const Reconstruction& p){
    sizeX =p.sizeX;
    sizeY =p.sizeY;
    value.resize(p.value.size());
    for(int i=0; i<p.value.size(); i++){
      value[i]=p.value[i];
    }
  }

  Reconstruction(Reconstruction&& p){
    sizeX =p.sizeX;
    sizeY =p.sizeY;
    value =std::move(p.value);
  }


  Reconstruction& operator=(const Reconstruction& p) {
    sizeX =p.sizeX;
    sizeY =p.sizeY;
    value.resize(p.value.size());
    for(int i=0; i<p.value.size(); i++){
      value[i]=p.value[i];
    }
    return *this;
  }

  Reconstruction& operator=(Reconstruction&& p){
    sizeX =p.sizeX;
    sizeY =p.sizeY;
    value =std::move(p.value);
    return *this;
  }

  Reconstruction scale(const Real_t scalar) const {
    Reconstruction ret=*this;
    for(auto &v : ret.value ){
      v *=scalar;
    }
    return std::move(ret);
  }

  Reconstruction subtract(const Reconstruction &r) const {
    Reconstruction ret=*this;
    assert(ret.value.size()==r.value.size());
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]-=r.value[i];
    }
    return std::move(ret);
  }

  Reconstruction add(const Reconstruction &r) const {
    Reconstruction ret=*this;
    assert(ret.value.size()==r.value.size());
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]+=r.value[i];
    }
    return std::move(ret);
  }

  Real_t operator[](int i){
    return value[i];
  }

  Real_t getValue(int x, int y) const{
    return value[x +sizeX*y];
  }

  void setValue(int x, int y, Real_t val)
  {
    value[x +sizeX*y] =val;
  }

  void clear(double val =0.0){
    for(auto &v : value ){
      v=val;
    }
  }

  void stat(){
    Real_t mean=0.0;
    Real_t max=-9999999;
    Real_t min=9999999;
    for(auto v : value){
      mean+=v;
      max =std::max(v,max);
      min =std::min(v,min);
    }
    mean /=(Real_t)value.size();
    std::cout << "mean:" << mean << " max:" << max << " min:" << min << std::endl;
  }
};

class ProjectionOperator
{
private:
  Geometry geo;

public:
  ProjectionOperator(const Geometry &geoIn) : geo(geoIn) {}
  ~ProjectionOperator(){}

  Projection forwardProjection(const Reconstruction &rec){
    Projection ret(geo);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int p=0; p< geo.projectionCount; p++){
      Real_t radian = (Real_t)p / (Real_t)geo.projectionCount * M_PI * 2;
      for(int x=0; x< geo.detectorPixelNum; x++){
	FP gen(geo, x, cosf(radian), sinf(radian));
	MatrixElement_t mat;
	Real_t           sum = 0.0F;
	while (1) {
	  if (!gen.GetMatrixElement(mat)) break;
	  sum += mat.length * rec.getValue(mat.index.x, mat.index.y);
	}
	ret.setValue(x, p, sum);
      }
    }
    return std::move(ret);
  }


  Reconstruction backwardProjection(const Projection &proj){
    Reconstruction ret(geo);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int p=0; p< geo.projectionCount; p++){
      Real_t radian = (Real_t)p / (Real_t)geo.projectionCount * M_PI * 2;
      for(int x=0; x< geo.detectorPixelNum; x++){
	FP gen(geo, x, cosf(radian), sinf(radian));
	MatrixElement_t mat;
	while (1) {
	  if (!gen.GetMatrixElement(mat)) break;
	  Real_t sum =ret.getValue(mat.index.x, mat.index.y);
	  sum += mat.length * proj.getValue(x, p);
	  ret.setValue(mat.index.x, mat.index.y, sum);
	}
      }
    }
    return std::move(ret);
  }

};
