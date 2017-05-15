#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################

CPPSRC += main.cpp
CPPSRC += cnpy.cpp
CPPSRC += libarg.cpp

CCOPT += -std=c++11 -g3 -fopenmp -lpthread -O1 -lz
LDOPT  = -pthread -lgomp -lz

TARGET = tst

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
INC_DIR += -I../ImageUtil
INC_DIR += -I../lib
INC_DIR += -I../cnpy

vpath %.cpp ../lib
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../DataGenerator
vpath %.cpp ../cnpy

vpath %.c ../lib

include ../lib/default_host.mk

run:
	./tst

clean-image:
	rm ./output/*
