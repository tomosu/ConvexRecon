#pragma once
#include "lodepng.h"
#include "ProjectionOperator.hpp"
#include "Geometry.hpp"

#include <vector>
#include <algorithm>

std::vector<Reconstruction> Png2Reconstruction(const char* filename, const Geometry &geo);

void Reconstruction2Png(const char* filename, const Geometry &geo, std::vector<Reconstruction> in);

