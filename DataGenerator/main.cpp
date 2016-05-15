#include "ProjectionGenerator.h"
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


/****************************************
 * main
 */
int main(int argc, char **argv)
{
    /******************************************
     * read argv
     */
    char   phantomPath[MBF]            = "../config/phantom.ini";
    char   geometryPath[MBF]           = "../config/geometry.ini";

    Real_t noiseIntensity              = 0.05;
    Real_t digitized_air               = 16383;
    int    isTransmission              = 1;

    Geometry geo(geometryPath);
    geo.Print();

    /****************************************
     * generate projection
     */
    Projection data =CalcProjectionValue(geo, phantomPath);

    /****************************************
     * add noise
     */
    AddNoise(geo, noiseIntensity, digitized_air, &data);
}
