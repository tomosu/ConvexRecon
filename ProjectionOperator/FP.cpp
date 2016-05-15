#include "FP.hpp"
#include "ProjectionOperator.hpp"

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>


FP::FP(const Geometry &inGeometry,
       int                x,
       Real_t             vCos,
       Real_t             vSin)

    : m_geometry(inGeometry)
{
    this->cosTheta  = vCos;
    this->sinTheta  = vSin;

    this->rayOrigin = LeftRotateAroundZ(m_geometry.getSourceOrigin(x), m_geometry.gantryCenter);
    this->detectorCenter   = LeftRotateAroundZ(m_geometry.getDetectorPos(x), m_geometry.gantryCenter);

    this->volumePlaneMinX = m_geometry.volumeCenter.x - m_geometry.volumePitch.x * (m_geometry.volumeSize.x / 2);
    this->volumePlaneMaxX = m_geometry.volumeCenter.x + m_geometry.volumePitch.x * (m_geometry.volumeSize.x / 2);
    this->volumePlaneMinY = m_geometry.volumeCenter.y - m_geometry.volumePitch.y * (m_geometry.volumeSize.y / 2);
    this->volumePlaneMaxY = m_geometry.volumeCenter.y + m_geometry.volumePitch.y * (m_geometry.volumeSize.y / 2);

    GetNormalizedRayDirection(
        this->elementRayDirection,
        this->detectorCenter);

    this->planeOfIncidence = GetPlaneOfIncidence(
        this->elementRayDirection);

    GetFirstIntersection(
        this->firstIntersection,
        this->planeOfIncidence,
        this->elementRayDirection);

    GetFirstIntersectVoxelIdx(
        this->currentIdx,
        this->planeOfIncidence,
        this->firstIntersection);

    SetIncrFlag(
        this->incr,
        this->elementRayDirection);

    InitCurrentSegment(
        this->currentSegment,
        this->elementRayDirection,
        this->incr,
        this->firstIntersection,
        this->planeOfIncidence);

    CalcNormalizedVoxelPitch(
        this->normalizedVoxelPitch,
        this->elementRayDirection);


}

Real2_t FP::SwapXY(const Real2_t &v)
{
    Real2_t ret;
    ret.x = v.y;
    ret.y = v.x;
    return ret;
}

Real2_t FP::LeftRotateAroundZ(const Real2_t &vin, const Real2_t &rotCenter)
{
  Real2_t v=vin -rotCenter;
    Real2_t vrot;
    vrot.x = v.x * this->cosTheta - v.y * this->sinTheta;
    vrot.y = v.x * this->sinTheta + v.y * this->cosTheta;
    return vrot +rotCenter;
}

void FP::GetNormalizedRayDirection(Real2_t       &out,
                                    const Real2_t &rayTarget)
{
    Real_t dirX = rayTarget.x - this->rayOrigin.x;
    Real_t dirY = rayTarget.y - this->rayOrigin.y;

    Real_t length = sqrtf(dirX * dirX + dirY * dirY);

    out.x = dirX / length;
    out.y = dirY / length;
}

void FP::SetIncrFlag(Int2_t         &incr,
                      const Real2_t &rayDirection)
{
    if (rayDirection.x >= 0) {
        incr.x = 1;
    } else {
        incr.x = -1;
    }

    if (rayDirection.y >= 0) {
        incr.y = 1;
    } else {
        incr.y = -1;
    }
}

int FP::GetPlaneOfIncidence(const Real2_t &rayDirection)
{
    int   plane = PLANE_NONE;
    Real_t t     = 0.0F;

    ////////////////////////////////////////////////////////////
    // Check X-Plane (LEFT)
    if (this->rayOrigin.x > this->volumePlaneMaxX) {
        if (std::abs(rayDirection.x) > 0) {
            t = std::abs((this->volumePlaneMaxX - this->rayOrigin.x) / rayDirection.x);
        } else {
            t = FLT_MAX;
        }

        Real_t y = this->rayOrigin.y + rayDirection.y * t;
        if (  y <= this->volumePlaneMaxY
           && y >= this->volumePlaneMinY) {
            plane = PLANE_LEFT;
        }
    }

    ////////////////////////////////////////////////////////////
    // Check X-Plane (RIGHT)
    if (this->rayOrigin.x < this->volumePlaneMinX) {
        if (std::abs(rayDirection.x) > 0) {
            t = std::abs((this->volumePlaneMinX - this->rayOrigin.x) / rayDirection.x);
        } else {
            t = FLT_MAX;
        }

        Real_t y = this->rayOrigin.y + rayDirection.y * t;
        if (  y <= this->volumePlaneMaxY
           && y >= this->volumePlaneMinY) {
            plane = PLANE_RIGHT;
        }
    }

    ////////////////////////////////////////////////////////////
    // Check Y-Plane (BOTTOM)
    if (this->rayOrigin.y > this->volumePlaneMaxY) {
        if (std::abs(rayDirection.y) > 0) {
            t = std::abs((this->volumePlaneMaxY - this->rayOrigin.y) / rayDirection.y);
        } else {
            t = FLT_MAX;
        }

        Real_t x = this->rayOrigin.x + rayDirection.x * t;
        if (  x <= this->volumePlaneMaxX
           && x >= this->volumePlaneMinX) {
            plane = PLANE_BOTTOM;
        }
    }

    ////////////////////////////////////////////////////////////
    // Check Y-Plane (TOP)
    if (this->rayOrigin.y < this->volumePlaneMinY) {
        if (std::abs(rayDirection.y) > 0) {
            t = std::abs((this->volumePlaneMinY - this->rayOrigin.y) / rayDirection.y);
        } else {
            t = FLT_MAX;
        }

        Real_t x = this->rayOrigin.x + rayDirection.x * t;
        if (  x <= this->volumePlaneMaxX
           && x >= this->volumePlaneMinX) {
            plane = PLANE_TOP;
        }
    }
    return plane;
}

void FP::GetFirstIntersection(Real2_t       &v,
                               const int       planeOfIncidence,
                               const Real2_t &rayDirection)
{
    v.x = 0;
    v.y = 0;
    Real_t t = 0;

    switch (planeOfIncidence) {
        case PLANE_TOP:
            if (std::abs(rayDirection.y) > 0) {
                t = std::abs((rayOrigin.y - this->volumePlaneMinY) / rayDirection.y);
            } else {
                t = FLT_MAX;
            }
            v.x = rayOrigin.x + rayDirection.x * t;
            v.y = this->volumePlaneMinY;
            break;

        case PLANE_BOTTOM:
            if (std::abs(rayDirection.y) > 0) {
                t = std::abs((rayOrigin.y - this->volumePlaneMaxY) / rayDirection.y);
            } else {
                t = FLT_MAX;
            }
            v.x = rayOrigin.x + rayDirection.x * t;
            v.y = this->volumePlaneMaxY;
            break;

        case PLANE_RIGHT:
            if (std::abs(rayDirection.x) > 0) {
                t = std::abs((rayOrigin.x - this->volumePlaneMinX) / rayDirection.x);
            } else {
                t = FLT_MAX;
            }
            v.x = this->volumePlaneMinX;
            v.y = rayOrigin.y + rayDirection.y * t;
            break;

        case PLANE_LEFT:
            if (std::abs(rayDirection.x) > 0) {
                t = std::abs((rayOrigin.x - this->volumePlaneMaxX) / rayDirection.x);
            } else {
                t = FLT_MAX;
            }
            v.x = this->volumePlaneMaxX;
            v.y = rayOrigin.y + rayDirection.y * t;
            break;

        default:
            break;
    }
}

void FP::GetFirstIntersectVoxelIdx(Int2_t         &v,
                                    const int       planeOfIncidence,
                                    const Real2_t &firstIntersection)
{
    v.x = 0;
    v.y = 0;

    switch (planeOfIncidence) {
        case PLANE_TOP:
            v.y = 0;
            v.x = (int)floor((firstIntersection.x - this->volumePlaneMinX) / m_geometry.volumePitch.x);
            break;

        case PLANE_BOTTOM:
            v.y = m_geometry.volumeSize.y - 1;
            v.x = (int)floor((firstIntersection.x - this->volumePlaneMinX) / m_geometry.volumePitch.x);
            break;

        case PLANE_RIGHT:
            v.x = 0;
            v.y = (int)floor((firstIntersection.y - this->volumePlaneMinY) / m_geometry.volumePitch.y);
            break;

        case PLANE_LEFT:
            v.x = m_geometry.volumeSize.x - 1;
            v.y = (int)floor((firstIntersection.y - this->volumePlaneMinY) / m_geometry.volumePitch.y);
            break;
    }
}

void FP::CalcNormalizedVoxelPitch(Real2_t       &pitch,
                                   const Real2_t &rayDirection)
{
    if (std::abs(rayDirection.x) > 0) {
        pitch.x = m_geometry.volumePitch.x / std::abs(rayDirection.x);
    } else {
        pitch.x = FLT_MAX;
    }

    if (std::abs(rayDirection.y) > 0) {
        pitch.y = m_geometry.volumePitch.y / std::abs(rayDirection.y);
    } else {
        pitch.y = FLT_MAX;
    }
}

void FP::InitCurrentSegment(Real2_t       &currentSegment,
                             const Real2_t &rayDirection,
                             const Int2_t   &incr,
                             const Real2_t &firstIntersection,
                             const int       planeOfIncidence)
{
    switch (planeOfIncidence) {
        case PLANE_TOP:
        case PLANE_BOTTOM:
            currentSegment.y = 0.0F;

            if (incr.x > 0) {
                if (std::abs(rayDirection.x) > 0) {
                    currentSegment.x = fmodf(firstIntersection.x - this->volumePlaneMinX, m_geometry.volumePitch.x) / std::abs(rayDirection.x);
                } else {
                    currentSegment.x = FLT_MAX;
                }
            } else {
                if (std::abs(rayDirection.x) > 0) {
                    currentSegment.x = fmodf(this->volumePlaneMaxX - firstIntersection.x, m_geometry.volumePitch.x) / std::abs(rayDirection.x);
                } else {
                    currentSegment.x = FLT_MAX;
                }
            }
            break;

        case PLANE_RIGHT:
        case PLANE_LEFT:
            currentSegment.x = 0.0F;

            if (incr.y > 0) {
                if (std::abs(rayDirection.y) > 0) {
                    currentSegment.y = fmodf(firstIntersection.y - this->volumePlaneMinY, m_geometry.volumePitch.y) / std::abs(rayDirection.y);
                } else {
                    currentSegment.y = FLT_MAX;
                }
            } else {
                if (std::abs(rayDirection.y) > 0) {
                    currentSegment.y = fmodf(this->volumePlaneMaxY - firstIntersection.y, m_geometry.volumePitch.y) / std::abs(rayDirection.y);
                } else {
                    currentSegment.y = FLT_MAX;
                }
            }
            break;
    }
}

bool FP::GetMatrixElement(MatrixElement_t &mat)
{
    if (this->planeOfIncidence == PLANE_NONE) {
        return false;
    }

    if (  (this->currentIdx.x < 0)
       || (this->currentIdx.y < 0)
       || (this->currentIdx.x >= this->m_geometry.volumeSize.x)
       || (this->currentIdx.y >= this->m_geometry.volumeSize.y)) {
        return false;
    }

    mat.index.x = this->currentIdx.x;
    mat.index.y = this->currentIdx.y;

    Real2_t boundaryDistance;
    boundaryDistance.x = std::abs(this->normalizedVoxelPitch.x - this->currentSegment.x);
    boundaryDistance.y = std::abs(this->normalizedVoxelPitch.y - this->currentSegment.y);

    if (boundaryDistance.x >= boundaryDistance.y) {
        mat.length              = boundaryDistance.y;
        this->currentIdx.y     += this->incr.y;
        this->currentSegment.y  = 0.0F;
        this->currentSegment.x += boundaryDistance.y;
    } else {
        mat.length              = boundaryDistance.x;
        this->currentIdx.x     += this->incr.x;
        this->currentSegment.x  = 0.0F;
        this->currentSegment.y += boundaryDistance.x;
    }
    return true;
}

#ifdef _TEST_FP

#include <gtest/gtest.h>

// SimpleTest/*{{{*/

class SimpleTest : public::testing::Test
{
public:
    Geometry geo;
    Real2_t   rayOrigin;
    Image_t    rec;
    Real_t      theta;

    virtual void SetUp()
    {
        geo.sod = 2.0F;
        geo.sid = 4.0F;

        geo.detector.size   = 10;
        geo.detectorCenterPos.x = 5.0F;
        geo.detectorPitch  = 1.0F;

        rayOrigin.x = 0;
        rayOrigin.y = -geo.sod;

        rec.size.x = 2;
        rec.size.y = 2;

        rec.center.x = 0.0F;
        rec.center.y = 0.0F;

        rec.pitch.x = 1.0F;
        rec.pitch.y = 1.0F;

        theta = 0.0F * M_PI / 180.0F;
    }

    virtual void TearDown()
    {
    }
};

TEST_F(SimpleTest, Test1) {
    FP gen(rayOrigin, geo, 4, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 1);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 1);
        if (mat.index.x == 0 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 1 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test2) {
    theta = 90.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 4, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 1);
        if (mat.index.x == 1 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 1 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test3) {
    theta = 180.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 4, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 0);
        if (mat.index.x == 1 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 0 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test4) {
    theta = 270.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 4, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 1);
        if (mat.index.x == 0 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 0 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test5) {
    FP gen(rayOrigin, geo, 5, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 0);
        if (mat.index.x == 0 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 1 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test6) {
    theta = 90.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 5, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 1);
        if (mat.index.x == 0 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 0 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test7) {
    theta = 180.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 5, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 1);
        ASSERT_FALSE(mat.index.x == 1 && mat.index.y == 1);
        if (mat.index.x == 0 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 1 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

TEST_F(SimpleTest, Test8) {
    theta = 270.0F * M_PI / 180.0F;
    FP gen(rayOrigin, geo, 5, rec, cosf(theta), sinf(theta));

    MatrixElement_t mat;
    while (1) {
        if (!gen.GetMatrixElement(mat)) break;
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 0);
        ASSERT_FALSE(mat.index.x == 0 && mat.index.y == 1);
        if (mat.index.x == 1 && mat.index.y == 0) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
        if (mat.index.x == 1 && mat.index.y == 1) ASSERT_TRUE(std::abs(sqrt(1.0 * 1.0 + (0.5 / 4) * (0.5 / 4))) == mat.length);
    }
}

/*}}}*/

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif // ifdef _TEST_BP
