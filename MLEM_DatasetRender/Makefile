#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################

CPPSRC += main.cpp
CPPSRC += SetupPhantom.cpp
CPPSRC += libini.cpp
CPPSRC += libarg.cpp
CPPSRC += Generator.cpp
CPPSRC += ProjectionGenerator.cpp
CPPSRC += FP.cpp

CCOPT += -std=c++11 -g3 -fopenmp -lpthread -O1
LDOPT  = -pthread -lgomp

TARGET = tst

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
INC_DIR += -I../ImageUtil
INC_DIR += -I../lib

vpath %.cpp ../lib
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../DataGenerator

vpath %.c ../lib

include ../lib/default_host.mk

run:
	./tst

clean-image:
	rm ./output/*
