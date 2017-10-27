#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################

CPPSRC += main.cpp
CPPSRC += lodepng.cpp
CPPSRC += cnpy.cpp
CPPSRC += libini.cpp
CPPSRC += PngReconConverter.cpp
CPPSRC += NpyReconConverter.cpp
CPPSRC += FP.cpp

#CSRC += mt19937ar.c

CCOPT += -std=c++11 -g3 -fopenmp -lpthread -O1
TARGET = tst

#LDOPT   += -lgomp

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../cnpy
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
INC_DIR += -I../lib

vpath %.cpp ../lib
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../cnpy
vpath %.cpp ../DataGenerator

vpath %.c ../lib

include ../lib/default_host.mk

run:
	./tst

clean-image:
	rm ./output/*
