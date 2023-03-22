#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <complex.h> 
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <time.h>  
#include <mpi.h>
#include"fftw_f.h"
#include"NSE.h"
#include"file_management.h"

#ifndef MPI_DEBUG_H
#define MPI_DEBUG_H
#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)
#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)
#endif





#ifndef VAR_DEFS          
#define VAR_DEFS 1


#ifndef VAR_DECLS
# define _DECL extern
# define _INIT(variable)
#else
# define _DECL
# define _INIT(variable)  = variable
#endif

_DECL  double PI _INIT(3.14159265358979323846264338327950288419716939937510582097494459);
#endif


extern int Xpt;
extern int Ypt;
extern int Zpt;



extern int N;
extern int n;
extern int Nsor;

extern int rank;
extern int size;

extern int ib;
extern int ie;


extern int File_number;

extern long double dx;
extern long double dy;
extern long double dz;
extern long double dt;

extern long double Lx;
extern long double Ly;
extern long double Lz;

extern long double LN;


extern long double w;


extern long double Reynoldsnumber;
extern long double Froudenumber  ;
extern long double Prandtlnumber ;




#define Unp1(i,j,k) (unp1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Un(i,j,k)   (un  [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Unm1(i,j,k) (unm1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Vnp1(i,j,k) (vnp1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vn(i,j,k)   (vn  [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vnm1(i,j,k) (vnm1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Wnp1(i,j,k) (wnp1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Wn(i,j,k)   (wn  [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Wnm1(i,j,k) (wnm1[Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Vtx(i,j,k) (vtx [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vty(i,j,k) (vty [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vtz(i,j,k) (vtz [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Vttx(i,j,k) (vttx [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vtty(i,j,k) (vtty [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vttz(i,j,k) (vttz [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Udx(i,j,k) (udx [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Udy(i,j,k) (udy [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Udz(i,j,k) (udz [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Vdx(i,j,k) (vdx [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vdy(i,j,k) (vdy [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Vdz(i,j,k) (vdz [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define NUn(i,j,k) (Nun [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define NVn(i,j,k) (Nvn [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define NWn(i,j,k) (Nwn [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define NU1(i,j,k) (Nu1 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define NV1(i,j,k) (Nv1 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define NW1(i,j,k) (Nw1 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Bm0(i,j,k) (bm0 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Bm1(i,j,k) (bm1 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Bm2(i,j,k) (bm2 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define Bp(i,j,k) (bp [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define WB0(i,j,k) (Wb0 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define WB1(i,j,k) (Wb1 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define WB2(i,j,k) (Wb2 [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define P(i,j,k)  (p [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])
#define Pn(i,j,k) (pn [Ypt*Zpt*(i-rank*(Xpt/size)) + Zpt*(j) + k])

#define A(i,k)  (a [(Zpt)*(i)  + k])

#define Sx(i,k)  (sx [(Zpt)*(i)  + k])
#define Sz(i,k)  (sz [(Zpt)*(i)  + k])


