#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################

CPPSRC += main.cpp
CPPSRC += SetupPhantom.cpp
CPPSRC += libini.cpp
CPPSRC += Generator.cpp
CPPSRC += ProjectionGenerator.cpp
CPPSRC += FP.cpp

#CSRC += mt19937ar.c

CCOPT += -std=c++11 -stdlib=libc++
TARGET = tst

#LDOPT   += -lgomp

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
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
