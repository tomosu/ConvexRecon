#include "FFT.h"
#include "RealTypes.hpp"


// utils
void FFT_2d(int sizeX, int sizeY,
	    std::vector< std::complex<Real_t> > signal_in,
	    std::vector< std::complex<Real_t> > *signal_out){

  std::vector< std::complex<Real_t> > out_tmp(sizeX*sizeY);

  // x direction fft
  for(int y=0; y<sizeY; y++){
    std::vector< std::complex<Real_t> > x_tmp_in(sizeX);
    std::vector< std::complex<Real_t> > x_tmp_out(sizeX);
    for(int x=0; x<sizeX; x++){
      x_tmp_in[x] =signal_in[x +sizeX*y];
    }
    FFT_1d(x_tmp_in, &x_tmp_out);
    for(int x=0; x<sizeX; x++){
      out_tmp[x +y*sizeX] =x_tmp_out[x];
    }
  }

  // y direction fft
  for(int x=0; x<sizeX; x++){
    std::vector< std::complex<Real_t> > y_tmp_in(sizeY);
    std::vector< std::complex<Real_t> > y_tmp_out(sizeY);
    for(int y=0; y<sizeY; y++){
      y_tmp_in[y] =out_tmp[x+sizeX*y];
    }
    FFT_1d(y_tmp_in, &y_tmp_out);
    for(int y=0; y<sizeY; y++){
      out_tmp[x +y*sizeX] =y_tmp_out[y];
    }
  }

  auto& F = *signal_out;
  for(int i =0; i<signal_in.size(); i++){
    F[i] =out_tmp[i];
  }
}


void iFFT_2d(int sizeX, int sizeY,
	     std::vector< std::complex<Real_t> > signal_in,
	     std::vector< std::complex<Real_t> > *signal_out){

  std::vector< std::complex<Real_t> > out_tmp(sizeX*sizeY);

  // x direction fft
  for(int y=0; y<sizeY; y++){
    std::vector< std::complex<Real_t> > x_tmp_in(sizeX);
    std::vector< std::complex<Real_t> > x_tmp_out(sizeX);
    for(int x=0; x<sizeX; x++){
      x_tmp_in[x] =signal_in[x +sizeX*y];
    }
    iFFT_1d(x_tmp_in, &x_tmp_out);
    for(int x=0; x<sizeX; x++){
      out_tmp[x +y*sizeX] =x_tmp_out[x];
    }
  }

  // y direction fft
  for(int x=0; x<sizeX; x++){
    std::vector< std::complex<Real_t> > y_tmp_in(sizeY);
    std::vector< std::complex<Real_t> > y_tmp_out(sizeY);
    for(int y=0; y<sizeY; y++){
      y_tmp_in[y] =out_tmp[x+sizeX*y];
    }
    iFFT_1d(y_tmp_in, &y_tmp_out);
    for(int y=0; y<sizeY; y++){
      out_tmp[x +y*sizeX] =y_tmp_out[y];
    }
  }

  auto& F = *signal_out;
  for(int i =0; i<signal_in.size(); i++){
    F[i] =out_tmp[i];
  }

}


void FFT_1d(std::vector< std::complex<Real_t> > signal_in,
	    std::vector< std::complex<Real_t> > *signal_out){

  int n_level;
  std::vector< int > ids;
  {
    n_level = lc_fft_calc_ids( signal_in.size(), &ids );
  }
  lc_fft( signal_in, ids, n_level, signal_out, false );
}


void iFFT_1d(std::vector< std::complex<Real_t> > signal_in,
	     std::vector< std::complex<Real_t> > *signal_out){

  int n_level;
  std::vector< int > ids;
  {
    n_level = lc_fft_calc_ids( signal_in.size(), &ids );
  }

  lc_inverse_fft( signal_in, ids, n_level, signal_out );
}


// fft funcs
int lc_fft_calc_ids(const int N,
		    std::vector< int >* pids){

  int n_level;
  {
    auto& i = n_level;
    for( i=0; i<64; ++i )
      if( N>>i == 1) break;
  }
  std::vector< int >& ids = *pids;
  ids.reserve( N );
  {
    ids.push_back( 0 );
    ids.push_back( 1 );
    for( int i=0; i<n_level-1; ++i ){
      auto sz = ids.size();
      for_each( ids.begin(), ids.end(), [](int& x){ x*=2; } );
      ids.insert( ids.end(), ids.begin(), ids.end() );
      auto it = ids.begin();
      std::advance( it, sz );
      for_each( it, ids.end(), [](int&x){ x+=1; } );
    }// i
  }
  return n_level;
}


void lc_fft(const std::vector< std::complex<Real_t> >& a,
	    const std::vector< int >& ids,
	    const int n_level,
	    std::vector< std::complex< Real_t > >* pout,
	    bool is_inverse){

  auto N = a.size();
  auto& F = *pout;
  {
    F.resize( N );
    for( int i=0; i<N; ++i )
      F[ i ] = a[ids[i]];
  }

  unsigned int po2 = 1;
  for( int i_level=1; i_level<=n_level; ++i_level ){
    po2<<=1;
    const int po2m = po2>>1;
    auto w = exp( std::complex<Real_t>(.0, 2*M_PI/(Real_t)po2) );
    w = is_inverse ? conj(w): w;
    auto ws = std::complex<Real_t>(1,0);
    for( int k=0; k<po2m; ++k ){
      for( int j=0; j<N; j+=po2 ){
	auto pa = &F[j+k];
	auto pb = &F[j+k+po2m];
	auto wfb = ws**pb;
	*pb = *pa - wfb;
	*pa += wfb;
      }// j
      ws *= w;
    }// k
  }// i_level
}


void lc_inverse_fft( const std::vector< std::complex<Real_t> >& a,
		     const std::vector< int >& ids, const int n_level,
		     std::vector< std::complex< Real_t > >* pout ){

  lc_fft( a, ids, n_level, pout, true );
  auto N = a.size();
  std::for_each( pout->begin(), pout->end(),
		 [N](std::complex<Real_t>& val){val/=N;} );
}
