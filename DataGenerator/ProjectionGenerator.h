#pragma once
#include "Generator.h"
#include "Phantom.h"
#include "GeometricTools.h"
#include "ProjectionOperator.hpp"

#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define MBF 256

Projection CalcProjectionValue(const Geometry &geo,
			       char           *phantomPath);


void AddNoise(const Geometry &geo,
	      Real_t          noiseIntensity,
	      Real_t          digitized_air,
	      Projection           *proj,
	      unsigned int seed);
