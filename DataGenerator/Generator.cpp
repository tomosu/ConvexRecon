/************************************************************
 *
 * 2D Projection Generator
 *
 ************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "libarg.h"
#include "Generator.h"

Real_t InnerProduct(Real2_t v0, Real2_t v1)
{
    return v0.x * v1.x + v0.y * v1.y;
}

inline Real_t Det(Real_t a, Real_t b, Real_t c)
{
    Real_t d = b * b - 4 * a * c;
    return d;
}

inline Real_t QuadraFormula_1(Real_t a, Real_t b, Real_t c)
{
    Real_t s = (-b - sqrt(Det(a, b, c))) / (2 * a); // smaller root
    return s;
}

inline Real_t QuadraFormula_2(Real_t a, Real_t b, Real_t c)
{
    Real_t s = (-b + sqrt(Det(a, b, c))) / (2 * a); // bigger root
    return s;
}

Real_t GetNodeDistance(const Real2_t &origin,
                       const Real2_t &dist,
                       const Real2_t &circleCenter,
                       Real_t            rad)
{
    // unit ray direction vector
    Real2_t vd;
    vd.x = dist.x - origin.x;
    vd.y = dist.y - origin.y;

    Real_t length = sqrt(vd.x * vd.x + vd.y * vd.y);
    vd.x = vd.x / length;
    vd.y = vd.y / length;

    // vector pointing center of ball from ray origin
    Real2_t v;
    v.x = circleCenter.x - origin.x;
    v.y = circleCenter.y - origin.y;

    // determine coefficient
    Real_t a = 1;
    Real_t b = -2.0 * InnerProduct(v, vd);
    Real_t c = InnerProduct(v, v) - rad * rad;

    if (Det(a, b, c) > 0.0) {
        Real_t t1 = QuadraFormula_1(a, b, c);
        Real_t t2 = QuadraFormula_2(a, b, c);
        return fabs(t1 - t2);
    } else {
        return 0.0;
    }
}

#define INSIDE  1
#define OUTSIDE 0

int IsInside(Real2_t point,
             Real2_t circleCenter,
             Real_t     rad)
{
    Real2_t vdiff;
    vdiff.x = point.x - circleCenter.x;
    vdiff.y = point.y - circleCenter.y;
    if (InnerProduct(vdiff, vdiff) < rad * rad) {
        return INSIDE;
    }
    return OUTSIDE;
}

int save_pixel_bin(char   *filepath,
                   Real_t *prj_pixel,
                   int     size)
{
    int   ret = 1;
    int   i;
    FILE *fprj;

    if ((fprj = fopen(filepath, "w")) == NULL) {
        printf("result file open error\n");
        return 0;
    }

    float *array = (float *)malloc(sizeof(float) * size);

    for (int i = 0; i < size; i++) {
        *(array + i) = (float)*(prj_pixel + i);
    }

    fwrite(array, sizeof(float), size, fprj);

    fclose(fprj);

    free(array);

    return ret;
}


Real2_t RotateBodyAroundZero(Real_t     theta,
                                Real2_t v)
{
    Real2_t vrot;
    vrot.x = v.x * cos(theta) - v.y * sin(theta);
    vrot.y = v.x * sin(theta) + v.y * cos(theta);
    return vrot;
}


Real2_t RightRotateBodyAroundCenter(Real_t     theta,
				    const Real2_t &vin,
				    const Real2_t &rotationCenter)
{
    Real2_t v =vin-rotationCenter;

    Real2_t vrot;
    vrot.x = v.x * cos(-theta) - v.y * sin(-theta);
    vrot.y = v.x * sin(-theta) + v.y * cos(-theta);
    return vrot +rotationCenter;
}



Real2_t LeftRotateBodyAroundCenter(Real_t     theta,
				   const Real2_t &vin,
				   const Real2_t &rotationCenter)
{
    Real2_t v =vin-rotationCenter;

    Real2_t vrot;
    vrot.x = v.x * cos(theta) - v.y * sin(theta);
    vrot.y = v.x * sin(theta) + v.y * cos(theta);
    return vrot +rotationCenter;
}


float real1(void)
{
    return ((float)rand() + 1.0) / ((float)RAND_MAX + 2.0);
}
