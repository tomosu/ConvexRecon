#pragma once
#include "cnpy.h"
#include "ProjectionOperator.hpp"
#include "Geometry.hpp"

#include <vector>
#include <algorithm>


Reconstruction Npy2Reconstruction(const char* filename, const Geometry &geo);

void Reconstruction2Npy(const char* filename, Reconstruction in);
