#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################
CPPSRC += main.cpp
CPPSRC += PngReconConverter.cpp
CPPSRC += GradientFilter.cpp
CPPSRC += ChambolleProjectionOperator.cpp
CPPSRC += MultiChannelChambolleProjectionOperator.cpp
CPPSRC += lodepng.cpp
CPPSRC += FP.cpp
CPPSRC += libini.cpp


CCOPT += -std=c++11

TARGET = tst


INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../ImageUtil
INC_DIR += -I../ChambolleProjectionOperator
INC_DIR += -I../lib

vpath %.cpp ../lib
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../ImageUtil
vpath %.cpp ../ChambolleProjectionOperator

vpath %.c ../lib

include ../lib/default_host.mk


clean-image:
	rm ./output/*
