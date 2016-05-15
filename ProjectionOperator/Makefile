#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################
CPPSRC += main.cpp
#CPPSRC += libini.cpp
#CPPSRC += libarg.cpp
#CPPSRC += SetupPhantom.cpp
#CPPSRC += Generator.cpp
CPPSRC += FP.cpp


#CSRC += mt19937ar.c

#CCOPT += -fopenmp -O2
CCOPT += -std=c++11

TARGET = tst

#LDOPT   += -lgomp

INC_DIR += -I.
INC_DIR += -I../include
INC_DIR += -I../lib

vpath %.cpp ../lib
vpath %.c ../lib

include ../lib/default_host.mk

run:
	./tst -core 16 \
	-iteration 300 \
	-noise_intensity 1.0 \
	-digitized_air 16383 \
	-phantom_setup_file "phantoms/phantom.ini" \
	-geometry_setup_file "config.ini" \
	-output_dir "output" \
	-output_phantom_by_text "phantom.txt" \
	-output_recon_by_text "recon.txt" \
	-convert_float_to_unsigned_short 1

clean-image:
	rm ./output/*
