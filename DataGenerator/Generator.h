#pragma once

#include "RealTypes.hpp"


Real_t InnerProduct(Real2_t, Real2_t);

Real_t GetNodeDistance(const Real2_t &,
                       const Real2_t &,
                       const Real2_t &,
                       Real_t);

Real2_t RotateBodyAroundZero(Real_t     theta,
			     const Real2_t v);

Real2_t RightRotateBodyAroundCenter(Real_t     theta,
				    const Real2_t &v,
				    const Real2_t &rotCenter);


Real2_t LeftRotateBodyAroundCenter(Real_t     theta,
				   const Real2_t &v,
				   const Real2_t &rotCenter);
