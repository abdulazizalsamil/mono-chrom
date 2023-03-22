
#define VAR_DECLS
#include"include/main.h"





int Xpt=128;
int Ypt=1;
int Zpt=128;


int N    = 100000;

int n;

int Nsor=10;

long double dt=0.001;



long double Lx=20.0*3.14159265358979323846264338327950288419716939937510582097494459;
long double Ly=1.0;
long double Lz=4.0*3.14159265358979323846264338327950288419716939937510582097494459;

long double LN=(2.0/5.0)*3.14159265358979323846264338327950288419716939937510582097494459;




long double w=1.2;

long double Reynoldsnumber=  1000.0;
long double Froudenumber;
long double Prandtlnumber =  7.2;


long double dx;
long double dy;
long double dz;

int rank;
int size;
int ib;
int ie;

int File_number;


int main(int argc, char **argv){
n=0;

Froudenumber =sqrtl((2.0/5.0)*3.14159265358979323846264338327950288419716939937510582097494459);

int devN = N/100;


dx=( (long double)1.0)  / ( (long double)Xpt-1);
dy=( (long double)1.0)  / ( (long double)Ypt);
dz=( (long double)1.0)  / ( (long double)Zpt-1);

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
 ib = rank * (Xpt / size)+1;
 ie = (rank + 1) * (Xpt / size)+1;


complex long double * un   = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * vn   = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * wn   = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * p    = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * Nun = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Nvn = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Nwn = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * bm0 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * bm1 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * bm2 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));

complex long double * Wb0 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Wb1 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));
complex long double * Wb2 = (complex long double *)malloc((Xpt/size+2)*Ypt*Zpt*sizeof(complex long double));




if(rank==0){
time_t now;
time(&now);
printf("%s", ctime(&now));

File_number=0;
}
MPI_Bcast( &File_number, 1, MPI_INT, 0, MPI_COMM_WORLD); 



mat2fft(un);mat2fft(vn);mat2fft(wn);
mat2fft(p);
mat2fft(Nun);mat2fft(Nvn);mat2fft(Nwn);
mat2fft(bm0);mat2fft(bm1);mat2fft(bm2);
mat2fft(Wb0);mat2fft(Wb1);mat2fft(Wb2);



clock_t begin = clock();
clock_t begtt = clock();

for ( n = File_number*N+1; n <= N*(File_number+1); ++n){

NSE_solver(Nun,Nvn,Nwn,un,vn,wn,p,bm0,bm1,bm2,Wb0,Wb1,Wb2);




  if( n % devN==0 ){
    fft2mat(un);fft2mat(vn);fft2mat(wn);
    fft2mat(p);fft2mat(bm0);
    matrix_save(un,vn,wn,p,bm0,(n/devN));
    Wtime((n/devN), begtt);
    mat2fft(un);mat2fft(vn);mat2fft(wn);
    mat2fft(p);mat2fft(bm0);
    begtt = clock();
  }



}

  MPI_Barrier(MPI_COMM_WORLD);
fft2mat(un);fft2mat(vn);fft2mat(wn);
fft2mat(p);
fft2mat(Nun);fft2mat(Nvn);fft2mat(Nwn);
fft2mat(bm0);fft2mat(bm1);fft2mat(bm2);
fft2mat(Wb0);fft2mat(Wb1);fft2mat(Wb2);
  MPI_Barrier(MPI_COMM_WORLD);
NL_save(Nun,Nvn,Nwn,bm0 ,bm1,bm2,Wb0, Wb1,Wb2,(n/devN));



if(rank==0){

long double time_spent = 0.0;
clock_t end = clock();
time_spent += (long double)(end - begin) / CLOCKS_PER_SEC;

int hours  =(int)time_spent/3600;
time_spent =time_spent-hours*3600;
int minutes=(int)(int)time_spent/60;
time_spent =time_spent-minutes*60;
int seconds=(int)time_spent;

printf("The elapsed time is %d:%d:%d \n", hours,minutes,seconds);
}




MPI_Finalize();
free(un);free(vn);free(wn);
free(p);
free(Nun);free(Nvn);free(Nwn);
free(bm0);free(bm1);free(bm2);
free(Wb0);free(Wb1);free(Wb2);


  return 0;
}
