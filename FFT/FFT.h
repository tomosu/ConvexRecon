#pragma once

// "fft.h"
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include "RealTypes.hpp"

// utils
void FFT_2d(int sizeX, int sizeY,
	    std::vector< std::complex<Real_t> > signal_in,
	    std::vector< std::complex<Real_t> > *signal_out);

void iFFT_2d(int sizeX, int sizeY,
	     std::vector< std::complex<Real_t> > signal_in,
	     std::vector< std::complex<Real_t> > *signal_out);

void FFT_1d(std::vector< std::complex<Real_t> > signal_in,
	    std::vector< std::complex<Real_t> > *signal_out);

void iFFT_1d(std::vector< std::complex<Real_t> > signal_in,
	     std::vector< std::complex<Real_t> > *signal_out);

// fft header
int lc_fft_calc_ids( const int N, std::vector< int >* pids );

void lc_fft(const std::vector< std::complex<Real_t> >& a,
	    const std::vector< int >& ids, const int n_level,
	    std::vector< std::complex<Real_t> >* pout, bool is_inverse=0 );

void lc_inverse_fft(const std::vector< std::complex<Real_t> >& a,
		    const std::vector<int>& ids,const int n_level,
		    std::vector< std::complex<Real_t> >* pout );
