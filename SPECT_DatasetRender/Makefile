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
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter50"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter50"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter100"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter100"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter200"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter200"
	./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0006" \
	-emit_ratio 0.1 \
	-atten_ratio 0.0006 \
	-max_iteration 101
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter50"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter50"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter100"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter100"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter200"
	mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter200"
	./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0003" \
        -emit_ratio 0.1 \
        -atten_ratio 0.0003 \
	-max_iteration 101
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_with_init/iter10"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_with_init/iter20"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_with_init/iter30"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_with_init/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_with_init/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_with_init/iter10"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_with_init/iter20"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_with_init/iter30"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_with_init/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_with_init/iter100"
	# ./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0006" \
	# -emit_ratio 0.1 \
	# -atten_ratio 0.0006 \
	# -max_iteration 101 \
	# -with_init_dir "../dcnn_results_for_init/e0.1_a0.0006_spect_ct" \
	# -init_gain 1.0
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_with_init/iter10"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_with_init/iter20"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_with_init/iter30"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_with_init/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_with_init/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_with_init/iter10"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_with_init/iter20"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_with_init/iter30"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_with_init/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_with_init/iter100"
	# ./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0003" \
        # -emit_ratio 0.1 \
        # -atten_ratio 0.0003 \
	# -max_iteration 101 \
	# -with_init_dir "../dcnn_results_for_init/e0.1_a0.0003_spect_ct" \
	# -init_gain 1.0
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter10"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter20"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter30"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter10"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter20"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter30"
	#./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0006" \
	#-emit_ratio 0.1 \
	#-atten_ratio 0.0006
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter10"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter20"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter30"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter10"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter20"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter30"
	#./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0003" \
        #-emit_ratio 0.1 \
        #-atten_ratio 0.0003
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM_smooth/iter200"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM_smooth/iter200"
	# ./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0006" \
	# -emit_ratio 0.1 \
	# -atten_ratio 0.0006
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter50"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter100"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM_smooth/iter200"
	# mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM_smooth/iter200"
	# ./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0003" \
        # -emit_ratio 0.1 \
        # -atten_ratio 0.0003
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter50"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter50"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter100"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter100"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/MLEM/iter200"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0006/OSEM/iter200"
	#./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0006" \
	#-emit_ratio 0.1 \
	#-atten_ratio 0.0006
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter50"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter50"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter100"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter100"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/MLEM/iter200"
	#mkdir -p "test_img_iter_recon_results/spect_rec_e0.1_a0.0003/OSEM/iter200"
	#./tst -output_dir "test_img_iter_recon_results/spect_rec_e0.1_a0.0003" \
        #-emit_ratio 0.1 \
        #-atten_ratio 0.0003



clean-image:
	rm ./output/*
