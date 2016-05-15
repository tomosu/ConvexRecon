#pragma once

#include "RealTypes.hpp"
#include "Geometry.hpp"

#define PLANE_TOP    0
#define PLANE_LEFT   1
#define PLANE_RIGHT  2
#define PLANE_BOTTOM 3
#define PLANE_NONE   4

#define ORDER_YX 1
#define ORDER_XY 2

class FP{

public:

  FP (const Geometry &geometry,
      int              x,
      Real_t             vCos,
      Real_t             vSin);

  bool GetMatrixElement(MatrixElement_t &mat);

private:

  Real2_t SwapXY(const Real2_t &v);

  Real2_t LeftRotateAroundZ(const Real2_t &v, const Real2_t &rotationCenter);

  void GetNormalizedRayDirection(Real2_t       &out,
				 const Real2_t &rayTarget);

  int GetPlaneOfIncidence(const Real2_t &rayDirection);

  void GetFirstIntersection(Real2_t       &v,
			    const int       planeOfIncidence,
			    const Real2_t &rayDirection);

  void GetFirstIntersectVoxelIdx(Int2_t         &v,
				 const int       planeOfIncidence,
				 const Real2_t &firstIntersection);

  void InitCurrentSegment(Real2_t       &currentSegment,
			  const Real2_t &rayDirection,
			  const Int2_t   &incr,
			  const Real2_t &firstIntersection,
			  const int       planeOfIncidence);

  void SetIncrFlag(Int2_t         &incr,
		   const Real2_t &rayDirection);

  void CalcNormalizedVoxelPitch(Real2_t       &pitch,
				const Real2_t &rayDirection);

  const Geometry &m_geometry;

  Real_t sinTheta;
  Real_t cosTheta;

  Real2_t elementRayDirection;

  Int2_t   incr;
  Real2_t normalizedVoxelPitch;

  Real2_t currentSegment;
  Int2_t   currentIdx;

  Real2_t rayOrigin;
  Real2_t detectorCenter;
  Real2_t firstIntersection;
  int      planeOfIncidence;

  Real_t volumePlaneMinX;
  Real_t volumePlaneMaxX;
  Real_t volumePlaneMinY;
  Real_t volumePlaneMaxY;
};
