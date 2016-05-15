#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################
CPPSRC += main.cpp
CPPSRC += ChambolleProjectionOperator.cpp
CPPSRC += GradientFilter.cpp
CPPSRC += FP.cpp
CPPSRC += libini.cpp


CCOPT += -std=c++11

TARGET = tst

#LDOPT   += -lgomp

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../ChambolleProjectionOperator
INC_DIR += -I../lib

vpath %.cpp ../lib
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../ChambolleProjectionOperator

vpath %.c ../lib

include ../lib/default_host.mk


clean-image:
	rm ./output/*
