#pragma once
#include "RealTypes.hpp"
#include <cmath>

inline
Real_t CoordinateMethodTri(const Real2_t &p1,
                          const Real2_t &p2,
                          const Real2_t &p3)
{
    float s;
    s = std::abs(p1.x * p2.y - p2.x * p1.y +
              p2.x * p3.y - p3.x * p2.y +
              p3.x * p1.y - p1.x * p3.y) * (1.0F / 2.0F);
    return s;
}

inline
Real_t CoordinateMethodQuad(const Real2_t &p1,
                           const Real2_t &p2,
                           const Real2_t &p3,
                           const Real2_t &p4)
{
    float s;
    s = std::abs(p1.x * p2.y - p2.x * p1.y +
              p2.x * p3.y - p3.x * p2.y +
              p3.x * p4.y - p4.x * p3.y +
              p4.x * p1.y - p1.x * p4.y) * (1.0F / 2.0F);
    return s;
}

inline
Real3_t GetMassCenterTri(const Real3_t &p1,
                              const Real3_t &p2,
                              const Real3_t &p3)
{
    Real3_t p;
    p.x = (p1.x + p2.x + p3.x) * (1.0F / 3.0F);
    p.y = (p1.y + p2.y + p3.y) * (1.0F / 3.0F);
    p.z = (p1.z + p2.z + p3.z) * (1.0F / 3.0F);

    return p;
}

inline
Real3_t GetMassCenterQuad(const Real3_t &p1,
                               const Real3_t &p2,
                               const Real3_t &p3,
                               const Real3_t &p4)
{
    Real3_t p;
    p.x = (p1.x + p2.x + p3.x + p4.x) * (1.0F / 4.0F);
    p.y = (p1.y + p2.y + p3.y + p4.y) * (1.0F / 4.0F);
    p.z = (p1.z + p2.z + p3.z + p4.z) * (1.0F / 4.0F);

    return p;
}

inline
Real3_t GetMidPoint(const Real3_t &p1,
		    const Real3_t &p2)
{
    Real3_t p;
    p.x = (p1.x + p2.x) * (1.0F / 2.0F);
    p.y = (p1.y + p2.y) * (1.0F / 2.0F);
    p.z = (p1.z + p2.z) * (1.0F / 2.0F);

    return p;
}
