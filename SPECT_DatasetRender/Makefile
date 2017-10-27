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
CPPSRC += NpyReconConverter.cpp
CPPSRC += FP.cpp
CPPSRC += FFT.cpp
CPPSRC += cnpy.cpp

CCOPT += -std=c++11 -g3 -fopenmp -lpthread -O1
LDOPT  = -pthread -lgomp

TARGET = tst

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../lib
INC_DIR += -I../cnpy
INC_DIR += -I../FFT
INC_DIR += -I../ProjectionOperator
INC_DIR += -I../DataGenerator
INC_DIR += -I../ImageUtil

vpath %.cpp ../lib
vpath %.cpp ../FFT
vpath %.cpp ../ProjectionOperator
vpath %.cpp ../DataGenerator
vpath %.cpp ../ImageUtil
vpath %.cpp ../cnpy

vpath %.c ../lib

include ../lib/default_host.mk

run:
#	mkdir "spect_rec_e0.1_a0.0006"
#	./tst -output_dir "spect_rec_e0.1_a0.0006" \
#	-emit_ratio 0.1 \
#	-atten_ratio 0.0006
	mkdir "spect_rec_e0.1_a0.0003"
	./tst -output_dir "spect_rec_e0.1_a0.0003" \
        -emit_ratio 0.1 \
        -atten_ratio 0.0003
#	mkdir "spect_rec_e0.01_a0.0001"
#	./tst -output_dir "spect_rec_e0.01_a0.0001" \
#        -emit_ratio 0.01 \
#        -atten_ratio 0.0001
	mkdir "spect_rec_e0.01_a0.0003"
	./tst -output_dir "spect_rec_e0.01_a0.0003" \
        -emit_ratio 0.01 \
        -atten_ratio 0.0003



clean-image:
	rm ./output/*
