#include <math.h>
#include <stdio.h>
#include<stdlib.h>
#include <inttypes.h>
#include <assert.h>


complex long double U_top( int i);
complex long double Nb( int k);

void IC( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2);
void NSE_solver( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2);
void Non_linear( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *vtx,complex long double *vty,complex long double *vtz,complex long double *bp);

void P_solver( complex long double *p,complex long double *vtx,complex long double *vty ,complex long double *vtz );

void Diffusion( complex long double *p,complex long double *vtx,complex long double *vty ,complex long double *vtz,complex long double *un,complex long double *vn,complex long double *wn );
void P_boundary( complex long double *p );
void U_boundary( complex long double *vtx,complex long double *vty,complex long double *vtz);

complex long double D2x(complex long double *un,int i, int j, int k);
complex long double D2z(complex long double *un,int i, int j, int k);
complex long double Dx (complex long double *un,int i, int j, int k);
complex long double Dz (complex long double *un,int i, int j, int k);



void thomalg( complex long double *a, complex long double *b, int Nm);
void A_p( complex long double *a,int j);
void A_u( complex long double *a,int j);

