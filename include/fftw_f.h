

#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <complex.h> 
#include<stdlib.h>
#include <inttypes.h>
#include <assert.h>



double kappa(int j);


void complex2FFT( complex long double *U);
void FFT2complex( complex long double *U);
void FFT_test( );


void mat2fft( complex long double *un);
void fft2mat( complex long double *un);



void  fft(fftwl_complex *in, fftwl_complex *out);
void ifft(fftwl_complex *in, fftwl_complex *out);