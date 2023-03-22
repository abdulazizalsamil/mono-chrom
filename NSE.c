#include"include/main.h"




complex long double Nb( int k){

  long double z_0=1.0/2.0;
 return pow(1.0/coshl( (Lz/LN)*(k*dz-z_0) ),2);

//if(k==Zpt/2||k==Zpt/2-1||k==Zpt/2-2||k==Zpt/2+1||k==Zpt/2+2){return 1.0;}
//else{return 0.0;}
}


void IC( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2){

zerou(un);zerou(vn);zerou(wn);
zerou(p);
zerou(Nun);zerou(Nvn);zerou(Nwn);
zerou(bm0);zerou(bm1);zerou(bm2);
zerou(Wb0);zerou(Wb1);zerou(Wb2);


MPI_Barrier(MPI_COMM_WORLD);
if(rank==0){
 char name[50];
    FILE *stream;
    sprintf(name, "analysis/files/parameters.txt");
  stream= fopen(name,"w");
     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf %5.30Lf %5.30Lf\n", ( long double)Xpt,( long double)Ypt,(long double)Zpt, (long double)N,dt);
     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf %5.30Lf %5.30Lf\n", Lx,Ly,Lz,creall(Reynoldsnumber),creall(Froudenumber));
     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf %5.30Lf %5.30Lf\n", LN,Ly,Lz,creall(Reynoldsnumber),creall(Froudenumber));
fclose(stream);

}
MPI_Barrier(MPI_COMM_WORLD);

matrix_save(un,vn,wn,p,bm0,0);
}


void NSE_solver( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2){

complex long double * vtx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vty = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vtz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
zerou(vtx);zerou(vty);zerou(vtz);
complex long double * bp = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
zerou(bp);
complex long double err=0.0;
complex long double bxx;
complex long double ubx; 
complex long double Kb=1.0/(Reynoldsnumber*Prandtlnumber); 




  for (int i = ib; i < ie; ++i){for (int j = 0; j < Ypt; ++j){for (int k = 0; k < Zpt; ++k){
      Bp(i,j,k)=(1.0/Froudenumber)*(1.0/Froudenumber)*(  (23.0/12.0)*Bm0(i,j,k) -(16.0/12.0)*Bm1(i,j,k) +(5.0/12.0)*Bm2(i,j,k)  );
    }}}


MPI_Barrier(MPI_COMM_WORLD);
Non_linear(Nun,Nvn,Nwn,un,vn,wn,vtx,vty,vtz,bp);

MPI_Barrier(MPI_COMM_WORLD);
P_solver(p,vtx,vty,vtz);

MPI_Barrier(MPI_COMM_WORLD);
Diffusion(p,vtx,vty,vtz,un,vn,wn);

MPI_Barrier(MPI_COMM_WORLD);
equation(bm2,bm1);equation(bm1,bm0);

for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
      Bm0(i,j,k)=  dt*( (23.0/12.0)*WB0(i,j,k) -(16.0/12.0)*WB1(i,j,k) +(5.0/12.0)*WB2(i,j,k)  ) +Bm1(i,j,k)  ;
    }}}
equation(Wb2,Wb1);equation(Wb1,Wb0);
MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts( bm0);

for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
      ubx=0.0;//( Un(i,j,k)*( Bm0((i+1),j,k) - Bm0((i-1),j,k) )/(2.0*Lx*dx) -I*kappa(j)*Vn(i,j,k)*Bm0(i,j,k) + Wn(i,j,k)*( ( Bm0(i,j,(k+1)) - Bm0(i,j,(k-1)) )/(2.0*Lz*dz) ) ) ;
      bxx=(Bm0((i+1),j,k) -2.0*Bm0(i,j,k)+ Bm0((i-1),j,k) )/pow(Lx*dx,2)-kappa(j)*kappa(j)*Bm0(i,j,k)+( Bm0(i,j,(k+1))   - 2.0*Bm0(i,j,k) +Bm0(i,j,(k-1)) )/pow(Lz*dz,2) ;
      WB0(i,j,k)= (Nb(k)*Nb(k))*Wn(i,j,k)- ubx+Kb*bxx ;
    }}}


/*
int Ib=ib;
int Ie=ie;
  if (rank ==  0    ){Ib=ib+1;}
  if (rank == size-1){Ie=ie-1;}
MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(un);fill_ghosts(vn);fill_ghosts(wn);
for (int i = Ib; i < Ie; ++i){
  for (int j = 1; j < Ypt-1; ++j){
    for (int k = 1; k < Zpt-1; ++k){
       err=err+abs(Dx(un,i,j,k) + (-I*kappa(j))*Vn(i,j,k) +Dz(wn,i,j,k));
    }}}
if( creall(err)>1e-5 ){printf("Error\n");exit(0);}

*/


free(vtx);free(vty);free(vtz);free(bp);
}

void Non_linear( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *vtx,complex long double *vty,complex long double *vtz,complex long double *bp){

complex long double * udx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * udy = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * udz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * vdx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vdy = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vdz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * Nu1 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Nv1 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Nw1 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

zerou(udx);zerou(udy);zerou(udz);
zerou(vdx);zerou(vdy);zerou(vdz);

equation(Nu1,Nun);equation(Nv1,Nvn);equation(Nw1,Nwn);
zerou(Nun);zerou(Nvn);zerou(Nwn);

//Wave absorbing layer
//*
complex long double B_0=1.0/dt;
complex long double x_0=1.0;
complex long double L_B=Lz/2.0;

for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
        NUn(i,j,k)= NUn(i,j,k) - dt*B_0*exp(-pow( ( Lx*(i*dx-x_0)/L_B) ,2)  )*Un(i,j,k);
        NVn(i,j,k)= NVn(i,j,k) - dt*B_0*exp(-pow( ( Lx*(i*dx-x_0)/L_B) ,2)  )*Vn(i,j,k);
        NWn(i,j,k)= NWn(i,j,k) - dt*B_0*exp(-pow( ( Lx*(i*dx-x_0)/L_B) ,2)  )*Wn(i,j,k);   
    }}}
//*/



int Ib=ib;
int Ie=ie;
  if (rank ==  0    ){Ib=ib+1;}
  if (rank == size-1){Ie=ie-1;}


 /*omit nonlinear
for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

        Udx(i,j,k)=(-I*kappa(j))*Un(i,j,k);
        Udy(i,j,k)=(-I*kappa(j))*Vn(i,j,k);
        Udz(i,j,k)=(-I*kappa(j))*Wn(i,j,k);

    }}}
fft2mat(udx);fft2mat(udy);fft2mat(udz);
fft2mat(un);fft2mat(vn);fft2mat(wn);

MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(un);fill_ghosts(vn);fill_ghosts(wn);
for (int i = Ib; i < Ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 1; k < Zpt-1; ++k){

        NUn(i,j,k)=NUn(i,j,k)+0.5*( Un(i,j,k)*Dx(un,i,j,k)+ Vn(i,j,k)*Udx(i,j,k) +Wn(i,j,k)*Dz(un,i,j,k) );
        NVn(i,j,k)=NVn(i,j,k)+0.5*( Un(i,j,k)*Dx(vn,i,j,k)+ Vn(i,j,k)*Udy(i,j,k) +Wn(i,j,k)*Dz(vn,i,j,k) );
        NWn(i,j,k)=NWn(i,j,k)+0.5*( Un(i,j,k)*Dx(wn,i,j,k)+ Vn(i,j,k)*Udz(i,j,k) +Wn(i,j,k)*Dz(wn,i,j,k) );

    }}}

zerou(udx);zerou(udy);zerou(udz);
for (int i = ib; i < ie ; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){


        Udx(i,j,k)=Vn(i,j,k)*Un(i,j,k);
        Udy(i,j,k)=Vn(i,j,k)*Vn(i,j,k);
        Udz(i,j,k)=Vn(i,j,k)*Wn(i,j,k);

        Vdx(i,j,k)=Un(i,j,k)*Un(i,j,k);
        Vdz(i,j,k)=Un(i,j,k)*Wn(i,j,k);

        Vdy(i,j,k)=Wn(i,j,k)*Wn(i,j,k);
    }}}

MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(vdx);fill_ghosts(udx);fill_ghosts(vdz);
for (int i = Ib; i < Ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 1; k < Zpt-1; ++k){

        NUn(i,j,k)=NUn(i,j,k)+0.5*( Dx(vdx,i,j,k) + Dz(vdz,i,j,k) );
        NVn(i,j,k)=NVn(i,j,k)+0.5*( Dx(udx,i,j,k) + Dz(udz,i,j,k) );
        NWn(i,j,k)=NWn(i,j,k)+0.5*( Dx(vdz,i,j,k) + Dz(vdy,i,j,k) );

    }}}

mat2fft(Nun);mat2fft(Nvn);mat2fft(Nwn);
mat2fft(udx);mat2fft(udy);mat2fft(udz);


for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

        NUn(i,j,k)= NUn(i,j,k)  +0.5*(-I*kappa(j))*Udx(i,j,k);
        NVn(i,j,k)= NVn(i,j,k)  +0.5*(-I*kappa(j))*Udy(i,j,k);
        NWn(i,j,k)= NWn(i,j,k)  +0.5*(-I*kappa(j))*Udz(i,j,k);
    }}}

mat2fft(un);mat2fft(vn);mat2fft(wn);

/// omit nonlinear*/

for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

        Vtx(i,j,k)= Un(i,j,k) +0.5*dt*(3.0*NUn(i,j,k) - NU1(i,j,k));
        Vty(i,j,k)= Vn(i,j,k) +0.5*dt*(3.0*NVn(i,j,k) - NV1(i,j,k));  
        Vtz(i,j,k)= Wn(i,j,k) +0.5*dt*(3.0*NWn(i,j,k) - NW1(i,j,k)) -dt*Bp(i,j,k);
    }}}


fft2mat(vtx);fft2mat(vty);fft2mat(vtz);
U_boundary(vtx,vty,vtz);
mat2fft(vtx);mat2fft(vty);mat2fft(vtz);

free(udx);free(udy);free(udz);
free(vdx);free(vdy);free(vdz);
free(Nu1);free(Nv1);free(Nw1);

}



void P_solver( complex long double *p,complex long double *vtx,complex long double *vty ,complex long double *vtz ){
complex long double * pn    = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * a = (complex long double *)malloc((Zpt)*(Zpt)*sizeof(complex long double));
complex long double * b = (complex long double *)malloc((Zpt)*sizeof(complex long double));
for (int k = 0; k < Zpt; ++k){b[k]=0.0;for (int k2 = 0; k2 < Zpt; ++k2){A(k,k2)=0.0;}}
int k;

fft2mat(p);
P_boundary(p);
mat2fft(p);

equation(pn,p);
complex long double De;
int Ib=ib;
int Ie=ie;
  if (rank ==  0    ){Ib=ib+1;}
  if (rank == size-1){Ie=ie-1;}

for (int nsor = 0; nsor < Nsor; ++nsor){
MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(p);fill_ghosts(pn);fill_ghosts(vtx);
  for (int j = 0; j < Ypt; ++j){
    for (int i = Ib; i < Ie; ++i){
     for (k = 1; k < Zpt-1; ++k){
      De=Dx(vtx,i,j,k)+(-I*kappa(j))*Vty(i,j,k)+Dz(vtz,i,j,k);
      b[k]=(pow(Lz*dz,2)/dt)*( De )-pow(Lz*dz,2)*P((i+1),j,k)/pow(Lx*dx,2)-pow(Lz*dz,2)*Pn((i-1),j,k)/pow(Lx*dx,2);
      }

     A_p(a,j);
     thomalg(a,b,Zpt);

     for (k = 1; k < Zpt-1; ++k){Pn(i,j,k)=P(i,j,k)+w*( b[k]-P(i,j,k) );}

     }
equation(p,pn);  }



}

fft2mat(p);
P_boundary(p);
mat2fft(p);

free(a);free(b);free(pn);
}

void Diffusion(complex long double *p,complex long double *vtx,complex long double *vty ,complex long double *vtz,complex long double *un,complex long double *vn,complex long double *wn){

complex long double * udx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * udy = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * udz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * vdx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vdy = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vdz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * vttx = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vtty = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vttz = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * a = (complex long double *)malloc((Zpt)*(Zpt)*sizeof(complex long double));
complex long double * b = (complex long double *)malloc((Zpt)*sizeof(complex long double));

for (int k = 0; k < Zpt; ++k){b[k]=0.0;for (int k2 = 0; k2 < Zpt; ++k2){A(k,k2)=0.0;}}
int jb;int je;int k;int i;
complex long double Re=1.0/Reynoldsnumber;
int Ib=ib;
int Ie=ie;
  if (rank ==  0    ){Ib=ib+1;}
  if (rank == size-1){Ie=ie-1;}
zerou(vttx);zerou(vtty);zerou(vttz);
zerou(udx);zerou(udy);zerou(udz);
zerou(vdx);zerou(vdy);zerou(vdz);
MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(p);
for (int i = Ib; i < Ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (   k = 1; k < Zpt-1; ++k){

        Vttx(i,j,k)= -1.0*Dx(p,i,j,k)*dt+Vtx(i,j,k);
        Vtty(i,j,k)= -(-I*kappa(j) )*P(i,j,k)*dt+Vty(i,j,k);
        Vttz(i,j,k)= -1.0*Dz(p,i,j,k)*dt+Vtz(i,j,k);
    }}}

fft2mat(vtx);fft2mat(vty);fft2mat(vtz);
U_boundary(vtx,vty,vtz);
mat2fft(vtx);mat2fft(vty);mat2fft(vtz);


zerou(udx);zerou(udy);zerou(udz);


MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(un);fill_ghosts(vn);fill_ghosts(wn);
for (int i = Ib; i < Ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (   k = 1; k < Zpt-1; ++k){
        Vdx(i,j,k)=D2x(un,i,j,k)+(-kappa(j)*kappa(j))*Un(i,j,k)+D2z(un,i,j,k);
        Vdy(i,j,k)=D2x(vn,i,j,k)+(-kappa(j)*kappa(j))*Vn(i,j,k)+D2z(vn,i,j,k);
        Vdz(i,j,k)=D2x(wn,i,j,k)+(-kappa(j)*kappa(j))*Wn(i,j,k)+D2z(wn,i,j,k);

    }}}


fft2mat(vttx);fft2mat(vtty);fft2mat(vttz);
U_boundary(vttx,vtty,vttz);
mat2fft(vttx);mat2fft(vtty);mat2fft(vttz);

fft2mat(vdx);fft2mat(vdy);fft2mat(vdz);
U_boundary(vdx,vdy,vdz);
mat2fft(vdx);mat2fft(vdy);mat2fft(vdz);



for (int nsor = 0; nsor < Nsor; ++nsor){

MPI_Barrier(MPI_COMM_WORLD);
fill_ghosts(un);fill_ghosts(udx);
fill_ghosts(vn);fill_ghosts(udy);
fill_ghosts(wn);fill_ghosts(udz);
  for (int j = 0; j < Ypt; ++j){
    for ( i = Ib; i < Ie; ++i){
    



      A_u(a,j);
      for (k = 0; k < Zpt; ++k){
        b[k]=(0.5*dt*Re)*Vdx(i,j,k) + Vttx(i,j,k) + (0.5*dt*Re)*(  (Udx((i-1),j,k)+Un((i+1),j,k))/pow(Lx*dx,2)   );
      }
      thomalg(a,b,Zpt);
      for (k = 1; k < Zpt-1; ++k){Udx(i,j,k)=Un(i,j,k)+w*( b[k]-Un(i,j,k) );}

      A_u(a,j);
      for (k = 0; k < Zpt; ++k){
        b[k]=(0.5*dt*Re)*Vdy(i,j,k) + Vtty(i,j,k) + (0.5*dt*Re)*(  (Udy((i-1),j,k)+Vn((i+1),j,k))/pow(Lx*dx,2)   );
      }
      thomalg(a,b,Zpt);
      for (k = 1; k < Zpt-1; ++k){Udy(i,j,k)=Vn(i,j,k)+w*( b[k]-Vn(i,j,k) );}


      A_u(a,j);
      for (k = 0; k < Zpt; ++k){
        b[k]=(0.5*dt*Re)*Vdz(i,j,k) + Vttz(i,j,k) + (0.5*dt*Re)*(  (Udz((i-1),j,k)+Wn((i+1),j,k))/pow(Lx*dx,2)   );
      }
      thomalg(a,b,Zpt);
      for (k = 1; k < Zpt-1; ++k){Udz(i,j,k)=Wn(i,j,k)+w*( b[k]-Wn(i,j,k) );}

    }

equation(un,udx);equation(vn,udy);equation(wn,udz);

  }
}



fft2mat(un);fft2mat(vn);fft2mat(wn);
U_boundary(un,vn,wn);
mat2fft(un);mat2fft(vn);mat2fft(wn);

free(vttx);free(vtty);free(vttz);
free(udx);free(udy);free(udz);
free(vdx);free(vdy);free(vdz);
free(a);free(b);
}





void P_boundary( complex long double *p){



   if (rank == 0    ){
      for (int j = 0; j < Ypt; ++j){
      for (int  k = 1; k < Zpt-1; ++k){
         P(   1   ,j,k) = P(  2  ,j,k);
         P(   0   ,j,k) = P(  1  ,j,k);
        }}

  }
   if (rank == size-1){

      for (int j = 0; j < Ypt; ++j){
      for (int  k = 1; k < Zpt-1; ++k){
         P(( Xpt ),j,k) = P((Xpt-1),j,k);
         P((Xpt+1),j,k) = P((Xpt-1),j,k);
        }}
   }


for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){

    P(i,j,   0   ) = P( i,j,  1   );
    P(i,j,(Zpt-1)) = P(i,j,(Zpt-2));

    }}





}



void U_boundary( complex long double *vtx,complex long double *vty,complex long double *vtz){


complex long double * wm   = (complex long double *)malloc(Zpt*sizeof(complex long double));
complex long double * um   = (complex long double *)malloc(Zpt*sizeof(complex long double));
for (int k = 0; k < Zpt; ++k){wm[k]=0.0;um[k]=0.0;}


long double z_0     = 1.0/2.0;
long double A_tilda =0.1*(1.0-exp(-pow(dt*n/2.0,2) ) );
long double Bz,DBz;
    for (int k = 1; k < Zpt-1; ++k){
    
    Bz = A_tilda*exp(-  pow(( Lz*(k*dz-z_0)/LN),2)  );
    DBz= A_tilda*exp(-  pow(( Lz*(k*dz-z_0)/LN),2)  )*(Lz/(LN*LN))*(k*dz);


    wm[k]=Bz *sin( LN*dt*n );
    um[k]=DBz*cos( LN*dt*n );
}








/*
  

int         N1=32;
long double x_b=Lx/2.0;
long double z_0=Lz/2.0;
long double lampda_bar=0.9*(Lz/2.0);

long double sigma_bar=(LN/(Froudenumber*Froudenumber)) *(2*PI)/(lampda_bar);


long double delta_sigma=(2.0/3.0)*(sigma_bar/(N1-1.0));

long double k_bar=2.0*PI/lampda_bar;
long double cg=0.5*(sigma_bar)/(k_bar);


long double sigma_j=lampda_bar/20.0;

long double k_j=pow(sigma_j,2)/(LN/(Froudenumber*Froudenumber));


long double Bz,DBz;
long double Btime=(k_j/sigma_j);

long double A_tilda =0.1*(1.0-exp(-pow(dt*n/2.0,2) ) );


  for (int n1 = 0; n1 < N1; ++n1){
    sigma_j=sigma_j+delta_sigma; k_j=pow(sigma_j,2)/(LN/(Froudenumber*Froudenumber));
  }

long double Ntime=((k_j/sigma_j)-Btime)*2.0*x_b;




sigma_j=lampda_bar/20.0;
k_j=pow(sigma_j,2)/(LN/(Froudenumber*Froudenumber));
Btime=dt*n;


if(Btime<Ntime){
  for (int n1 = 0; n1 < N1; ++n1){
       sigma_j=sigma_j+delta_sigma; k_j=pow(sigma_j,2)/(LN/(Froudenumber*Froudenumber));
    for (int k = 1; k < Zpt-1; ++k){
    
    Bz = exp(-  pow(( (k*Lz*dz-z_0)/LN),2)  )* A_tilda ;
    DBz=-((k*Lz*dz-0.5)*(2.0/(LN*LN)))*exp(-  pow(( (k*Lz*dz-z_0)/LN),2)  )* A_tilda ;
    wm[k]=wm[k]+Bz*( sigma_j*( sin( (sigma_j/cg-k_j)*x_b-sigma_j*Btime ) )  );
    um[k]=um[k]+DBz*((sigma_j/k_j)*( cos( (sigma_j/cg-k_j)*x_b-sigma_j*Btime ) )  );
  }}
}

//*/


for (int i = ib-1; i < ie+1; ++i){
  for (int j = 0; j < Ypt; ++j){
        Vtx(i,j,(  0  ))= 0.0;Vty(i,j,(  0  ))= 0.0;Vtz(i,j,(  0  ))= 0.0;
        Vtx(i,j,(Zpt-1))= 0.0;Vty(i,j,(Zpt-1))= 0.0;Vtz(i,j,(Zpt-1))= 0.0;
    }}


   if (rank == 0    ){
      for (int j = 0; j < Ypt; ++j){
      for (int k = 0; k < Zpt; ++k){
        Vtx(1,j,k)=um[k];
        Vty(1,j,k)=0.0;
        Vtz(1,j,k)=wm[k];
        Vtx(0,j,k)=Vtx(( 1 ),j,k);
        Vty(0,j,k)=Vty(( 1 ),j,k);
        Vtz(0,j,k)=Vtz(( 1 ),j,k);
    }}
  }if (rank == size-1){

      for (int j = 0; j < Ypt; ++j){
      for (int k = 0; k < Zpt; ++k){
        Vtx(( Xpt ),j,k)=0.0;
        Vty(( Xpt ),j,k)=0.0;
        Vtz(( Xpt ),j,k)=0.0;
        Vtx((Xpt+1),j,k)=Vtx(( Xpt ),j,k);
        Vty((Xpt+1),j,k)=Vty(( Xpt ),j,k);
        Vtz((Xpt+1),j,k)=Vtz(( Xpt ),j,k);
    }}

   }




free(wm);free(um);

}


complex long double D2x(complex long double *un,int i, int j, int k){
return ( Un((i+1),j,k) - 2.0*Un(i,j,k) +Un((i-1),j,k) )/pow(Lx*dx,2);

}


complex long double D2z(complex long double *un,int i, int j, int k){
return ( Un(i,j,(k+1))   - 2.0*Un(i,j,k) +Un(i,j,(k-1)) )/pow(Lz*dz,2);

}

complex long double Dx(complex long double *un,int i, int j, int k){
return ( Un((i+1),j,k) - Un((i-1),j,k) )/(2.0*Lx*dx);


}

complex long double Dz(complex long double *un,int i, int j, int k){
return ( Un(i,j,(k+1)) - Un(i,j,(k-1)) )/(2.0*Lz*dz);

}

void thomalg( complex long double *a,complex long double *b, int Nm){
complex long double r=0.0;

int k;
 for ( k = 1; k < Nm; ++k){
  r=A((k-1),(k))/A((k-1),(k-1));
  A(k,k)=A((k),(k))-r*A((k),(k-1));
  b[k]  = b[k]-r*b[k-1];
  }
k=Nm-1;
  b[k]  = b[k]/A(k,k);


   for (int k = Nm-2; k >=0 ; --k){
    b[k]=(b[k]-A((k+1),k)*b[k+1])/A(k,k);
    }
}


void A_p( complex long double *a,int j){

 int k=1;
 A((k-1),(k-1))=-1.0*pow(Lz*dz,2)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2));A((k),(k-1))=1.0;
  for (k = 2; k < Zpt; ++k){A((k-1),(k-1))=-1.0*pow(Lz*dz,2)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2));A((k+1-1),(k-1))=1.0;A((k-1-1),(k-1))=1.0;}
 k=Zpt-1;
A(k,k)=-1.0*pow(Lz*dz,2)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2));A((k-1),k)=1.0;
}

void A_u( complex long double *a,int j){
  complex long double Re=1.0/Reynoldsnumber;
  int k=1;
  A((k-1),(k-1))=(1.0+(0.5*dt*Re)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2)));A((k),(k-1))=-(0.5*dt*Re)/pow(Lz*dz,2);
  for (k = 2; k < Zpt; ++k){
  A((k-1),(k-1))=(1.0+(0.5*dt*Re)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2)));A((k+1-1),(k-1))=-(0.5*dt*Re)/pow(Lz*dz,2);A((k-1-1),(k-1))=-(0.5*dt*Re)/pow(Lz*dz,2);}
  k=Zpt-1;
  A(k,k)=(1.0+(0.5*dt*Re)*(2.0/pow(Lx*dx,2)+kappa(j)*kappa(j)+2.0/pow(Lz*dz,2)));A((k-1),k)=-(0.5*dt*Re)/pow(Lz*dz,2);
}
