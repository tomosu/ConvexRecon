#pragma once
#include <cmath>

using Real_t = double;

struct Real2_t {
  Real2_t(Real_t xin=0, Real_t yin=0){x=xin, y=yin;}
  Real_t x;
  Real_t y;
  Real2_t operator+(const Real2_t &rhs) const { return Real2_t(x+rhs.x, y+rhs.y); }
  Real2_t operator-(const Real2_t &rhs) const { return Real2_t(x-rhs.x, y-rhs.y); }
  Real2_t scale(const Real_t scalar) const { return Real2_t(x*scalar, y*scalar);}
  Real_t norm_L1(){return std::abs(x) +std::abs(y);}
  Real_t norm_L2(){return sqrt(x*x +y*y);}

};


struct Real3_t {
  Real_t x;
  Real_t y;
  Real_t z;
};

struct Int2_t {
  int x;
  int y;
};

struct MatrixElement_t {
  Int2_t index;
  float length;
};
