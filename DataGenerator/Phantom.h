#pragma once
#include "RealTypes.hpp"

#define BASE     1
#define NOT_BASE 0

typedef struct
{
    Real_t     rad;
    Real_t     mu;
    Real2_t center;
    int        isBase;
} Phantom_t;

Phantom_t *setupPhantom(char *);

int count_body(char *);
