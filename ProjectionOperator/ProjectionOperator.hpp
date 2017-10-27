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


  Projection add_scalar(const Real_t v) const {
    Projection ret=*this;
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]+=v;
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


  Reconstruction add_scalar(const Real_t v) const {
    Reconstruction ret=*this;
    for(int i=0; i<ret.value.size(); i++){
      ret.value[i]+=v;
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



class EmissionProjectionOperator
{
private:
  Geometry geo;

public:

  struct RayPath{
    RayPath() :Attenuation(0.0), Emission(0.0), Length(0.0) {}
    RayPath(Real_t atten_in, Real_t emit_in, Real_t len_in)
      :Attenuation(atten_in), Emission(emit_in), Length(len_in)
    {}
    Real_t Attenuation;
    Real_t Emission;
    Real_t Length;
  };

  EmissionProjectionOperator(const Geometry &geoIn) : geo(geoIn) {}
  ~EmissionProjectionOperator(){}

  Real_t path2count(std::vector<RayPath> path, const Real_t emitRatio, const Real_t attenRatio){
    Real_t count =0.0;
    for(int i=0; i<path.size(); i++){
      Real_t emit =path[i].Emission *emitRatio;
      Real_t atten =0.0;
      for(int j=i; j<path.size(); j++){
	atten += path[j].Attenuation *path[j].Length;
      }
      count += emit *exp(-atten *attenRatio);
    }
    return count;
  }


  Projection forwardProjectionWithAttenuation(const Reconstruction &recEmit, const Reconstruction &recAtten,
					      const Real_t emitRatio, const Real_t attenRatio){
    Projection ret(geo);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int p=0; p< geo.projectionCount; p++){
      Real_t radian = (Real_t)p / (Real_t)geo.projectionCount * M_PI * 2;
      for(int x=0; x< geo.detectorPixelNum; x++){
	std::vector<RayPath> path;
	FP gen(geo, x, cosf(radian), sinf(radian));
	MatrixElement_t mat;
	Real_t          sum = 0.0F;
	while (1) {
	  if (!gen.GetMatrixElement(mat)) break;
	  path.push_back(RayPath(recAtten.getValue(mat.index.x, mat.index.y),
				 recEmit.getValue(mat.index.x, mat.index.y),
				 mat.length));
	}
	ret.setValue(x, p, path2count(path, emitRatio, attenRatio));
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
