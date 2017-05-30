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
CPPSRC += FFT.cpp

CCOPT += -std=c++11 -g3 -fopenmp -lpthread -O1
LDOPT  = -pthread -lgomp

TARGET = tst

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../lib
INC_DIR += -I../FFT
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
INC_DIR += -I../ImageUtil

vpath %.cpp ../lib
vpath %.cpp ../FFT
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../DataGenerator

vpath %.c ../lib

include ../lib/default_host.mk

run:
	./tst

prepare_dir:
	mkdir results
	mkdir results/answer
	mkdir results/output
	mkdir results/output/iter25
	mkdir results/output/iter50
	mkdir results/output/iter100
	mkdir results/output/iter150
	mkdir results/output/iter200
	mkdir results/output_fbp
	mkdir results/noise_output
	mkdir results/noise_output/iter25
	mkdir results/noise_output/iter50
	mkdir results/noise_output/iter100
	mkdir results/noise_output/iter150
	mkdir results/noise_output/iter200
	mkdir results/noise_output_fbp


clean-image:
	rm ./output/*
