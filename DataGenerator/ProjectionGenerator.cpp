#include "Generator.h"
#include "Phantom.h"
#include "GeometricTools.h"
#include "ProjectionOperator.hpp"

#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <random>

#define MBF 256


/****************************************
 * generate projection
 */
Projection CalcProjectionValue(const Geometry &geo,
			       char           *phantomPath)
{

  printf("read phantom ");

  int        phantomNum = count_body(phantomPath);
  Phantom_t *phantoms   = setupPhantom(phantomPath);

  printf("phantomNum: %d\n", phantomNum);
  printf("generate projection ");

  Projection ret(geo);
  Real_t omega = 2 * M_PI / geo.projectionCount;
  for (int prjID = 0; prjID < geo.projectionCount; prjID++) {
    for (int detID = 0; detID < geo.detectorPixelNum; detID++) {

      Real_t  theta           = omega * prjID;                       // rotate angle
      Real2_t rayOrigin =geo.getSourceOrigin(detID);
      Real2_t rayOriginRotate = LeftRotateBodyAroundCenter(theta, rayOrigin, geo.gantryCenter); // rotate
      Real2_t detectorCenter =geo.getDetectorPos(detID);
      Real2_t detectorCenterRotate = LeftRotateBodyAroundCenter(theta, detectorCenter, geo.gantryCenter);   // rotate

      Real_t atten_sum = 0.0;
      for (int i = 0; i < phantomNum; i++) {
	if (phantoms[i].isBase == NOT_BASE) {
	  // atten by signal
	  atten_sum += GetNodeDistance(rayOriginRotate,
				       detectorCenterRotate,
				       phantoms[i].center,
				       phantoms[i].rad) * phantoms[i].mu;
	} else {
	  // atten by base phatom
	  Real_t signal_length = 0.0;
	  for (int j = 0; j < phantomNum; j++) {
	    if (phantoms[j].isBase == NOT_BASE) {
	      signal_length += GetNodeDistance(rayOriginRotate,
					       detectorCenterRotate,
					       phantoms[j].center,
					       phantoms[j].rad);
	    }
	  }

	  Real_t base_length = GetNodeDistance(rayOriginRotate,
					       detectorCenterRotate,
					       phantoms[i].center,
					       phantoms[i].rad) - signal_length;

	  atten_sum += base_length * phantoms[i].mu;
	}
      }

      ret.setValue(detID, prjID, atten_sum);
    }
  }
  printf("done\n");
  free(phantoms);

  return std::move(ret);
}


/****************************************
 * add noise
 */


void AddNoise(const Geometry &geo,
	      Real_t          noiseIntensity,
	      Real_t          digitized_air,
	      Projection  *proj,
	      unsigned int seed)
{
  std::mt19937 mt(seed);
  std::uniform_real_distribution<double> genrand1(0.0,1.0);

    for (int prjID = 0; prjID < geo.projectionCount; prjID++) {
        for (int detID = 0; detID < geo.detectorPixelNum; detID++) {
            while (1) {
                Real_t air   = digitized_air;
                Real_t SND   = sqrtf(-2 * logf(genrand1(mt)) ) * cosf(2 * M_PI * genrand1(mt)); // BoxMuller
		Real_t projVal = proj->getValue(detID, prjID);
                Real_t count = air * exp(-projVal) + sqrt(air * exp(-projVal)) * SND * noiseIntensity;

                if (count > 0) {
		  proj->setValue(detID, prjID, -log(count / air));
		  break;
                } else {
		  continue;
                }
            }
        }
    }
}
