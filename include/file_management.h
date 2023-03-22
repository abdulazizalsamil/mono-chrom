#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include <inttypes.h>
#include <assert.h>


void equation( complex long double *vn, complex long double *un);
void zerou( complex long double *un);


void matrix_save(complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0,int n);
void NL_save(complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2,int ni);
void Open_files( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2, int n);


void fill_ghosts( complex long double *un);

void Wtime(int  Nf, clock_t begtt);


void Time_test( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2);



