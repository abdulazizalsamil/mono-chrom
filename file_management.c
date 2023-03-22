#include"include/main.h"

void equation( complex long double *vn, complex long double *un){

for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

      Vn(i,j,k)=Un(i,j,k);

    }}}


}


void zerou( complex long double *un){

for (int i = ib-1; i < ie+1; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

      Un(i,j,k)=0.0 + 0.0 * I;

    }}}


}

/*
matrix_save is used to save the date in files: the format is x,y,z,Unp1,Vpn1,Wpn1.
*/

void matrix_save(complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0,int n)
{
  char name[50];
  char namp[50];
    FILE *stream;
    FILE *fp;
  char c;
  sprintf(name, "analysis/files/V%dP%d.txt", n,rank);
  stream= fopen(name,"w");
 MPI_Barrier(MPI_COMM_WORLD);


for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){

     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf "   , creall(Un(i,j,k)),creall(Vn(i,j,k)),creall(Wn(i,j,k)));
     fprintf(stream, "%5.30Lf  %5.30Lf\n"   , creall(P(i,j,k)), creall(Bm0(i,j,k)) );
    }}}

    fclose(stream);
MPI_Barrier(MPI_COMM_WORLD);
 
if(rank==0){
    sprintf(name, "analysis/files/U%d.txt", n);
  stream= fopen(name,"w+");

for (int Pros = 0; Pros < size; ++Pros){
fseek(stream, 0,SEEK_END);

  sprintf(namp, "analysis/files/V%dP%d.txt", n,Pros);
  fp = fopen(namp,"r");

      c = fgetc(fp);
    while (c != EOF)
    {
        fputc(c, stream);
        c = fgetc(fp);
    }

  fclose(fp);
  remove(namp);
}


fclose(stream);
}

 MPI_Barrier(MPI_COMM_WORLD);


}



void NL_save(complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2,int ni)
{

  char name[50];
  char namp[50];
    FILE *stream;
    FILE *fp;
  char c;
  sprintf(name, "analysis/files/V%dP%d.txt", n,rank);
  stream= fopen(name,"w");
 MPI_Barrier(MPI_COMM_WORLD);


for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf "   , creall(NUn(i,j,k)),creall(NVn(i,j,k)),creall(NWn(i,j,k)));
     fprintf(stream, "%5.30Lf %5.30Lf "           , creall(Bm1(i,j,k)),creall(Bm2(i,j,k)));
     fprintf(stream, "%5.30Lf %5.30Lf %5.30Lf\n"  , creall(WB0(i,j,k)),creall(WB1(i,j,k)),creall(WB2(i,j,k)));
      }}}

    fclose(stream);
MPI_Barrier(MPI_COMM_WORLD);

if(rank==0){
    sprintf(name, "analysis/files/N%d.txt", ni-1);
  stream= fopen(name,"w+");

for (int Pros = 0; Pros < size; ++Pros){
fseek(stream, 0,SEEK_END);

  sprintf(namp, "analysis/files/V%dP%d.txt", n,Pros);
  fp = fopen(namp,"r");

      c = fgetc(fp);
    while (c != EOF)
    {
        fputc(c, stream);
        c = fgetc(fp);
    }

  fclose(fp);
  remove(namp);
}


fclose(stream);
}

 MPI_Barrier(MPI_COMM_WORLD);



  }



void Open_files( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2,int ni){

long double * Ub = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * Vb = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * Wb = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * Pb = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));

long double * B1 = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * B2 = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * W0 = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * W1 = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));
long double * W2 = (long double *)malloc(((Xpt+2)*Ypt*Zpt)*sizeof(long double));

zerou(un);zerou(vn);zerou(wn);
zerou(p);
zerou(Nun);zerou(Nvn);zerou(Nwn);
zerou(bm0);zerou(bm1);zerou(bm2);
zerou(Wb0);zerou(Wb1);zerou(Wb2);


for (int i = 0; i < Xpt+2; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
     Ub[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     Vb[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     Wb[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     Pb[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     B1[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     B2[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     W0[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     W1[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
     W2[Ypt*Zpt*(i) + Zpt*(j) + k]=0.0 ;
    }}}


  char name[45];
  FILE *stream;
  sprintf(name, "analysis/files/N%d.txt", ni);
  stream= fopen(name,"r");
  if (stream == NULL){printf("Error! opening file\n");exit(1);}

for (int i = 0; i < (Xpt*Ypt*Zpt); ++i){
       fscanf(stream, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf \n"   , &Ub[(i)],&Vb[(i)],&Wb[(i)], &B1[(i)],&B2[(i)],&W0[(i)],&W1[(i)],&W2[(i)]);
      }

  fclose(stream);


  for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
  for (int k = 0; k < Zpt; ++k){
      NUn((i),j,k)=Ub[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      NVn((i),j,k)=Vb[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      NWn((i),j,k)=Wb[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;

      Bm1((i),j,k)=B1[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      Bm2((i),j,k)=B2[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;

      WB0((i),j,k)=W0[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      WB1((i),j,k)=W1[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      WB2((i),j,k)=W2[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
    }}}


long double * Ua = (long double *)malloc((Xpt*Ypt*Zpt+2)*sizeof(long double));
long double * Va = (long double *)malloc((Xpt*Ypt*Zpt+2)*sizeof(long double));
long double * Wa = (long double *)malloc((Xpt*Ypt*Zpt+2)*sizeof(long double));
long double * Pa = (long double *)malloc((Xpt*Ypt*Zpt+2)*sizeof(long double));

long double * B0 = (long double *)malloc((Xpt*Ypt*Zpt+2)*sizeof(long double));


  sprintf(name, "analysis/files/U%d.txt", ni);
  stream= fopen(name,"r");
  if (stream == NULL){printf("Error! opening file\n");exit(1);}

for (int i = 0; i < (Xpt*Ypt*Zpt+2); ++i){

     fscanf(stream, "%Lf %Lf %Lf %Lf %Lf \n"   , &Ua[(i)],&Va[(i)],&Wa[(i)],&Pa[(i)],&B0[(i)]);
    }

    fclose(stream);


for (int i = ib; i < ie; ++i){
  for (int j = 0; j < Ypt; ++j){
    for (int k = 0; k < Zpt; ++k){
      Un ((i),j,k)=Ua[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      Vn ((i),j,k)=Va[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      Wn ((i),j,k)=Wa[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      P  ((i),j,k)=Pa[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
      Bm0((i),j,k)=B0[Ypt*Zpt*(i) + Zpt*(j) + k]+ 0.0 * I;
    }}}

free(Ua);free(Va);free(Wa);free(Pa);
free(Ub);free(Vb);free(Wb);free(Pb);
free(B1);free(B2);free(W0);free(W1);free(W2);
}


void fill_ghosts( complex long double *un){

 int tag_right_ghost = 100, tag_left_ghost = 101;


   int rank_left  = rank - 1;
   int rank_right = rank + 1;

   if (rank == 0    ){rank_left = size -1;}
   if (rank == size-1){ rank_right = 0;   }

   MPI_Request req[4];


   MPI_Isend(&Un( ib ,0,0)  , Ypt*Zpt, MPI_C_LONG_DOUBLE_COMPLEX, rank_left, tag_right_ghost , MPI_COMM_WORLD, &req[0]);
   MPI_Irecv(&Un(ib-1,0,0)  , Ypt*Zpt, MPI_C_LONG_DOUBLE_COMPLEX, rank_left, tag_left_ghost  , MPI_COMM_WORLD, &req[3]); 


   MPI_Irecv(&Un( ie ,0,0)  , Ypt*Zpt, MPI_C_LONG_DOUBLE_COMPLEX, rank_right, tag_right_ghost, MPI_COMM_WORLD, &req[1]);
   MPI_Isend(&Un(ie-1,0,0)  , Ypt*Zpt, MPI_C_LONG_DOUBLE_COMPLEX, rank_right, tag_left_ghost , MPI_COMM_WORLD, &req[2]);
   

   MPI_Waitall(4, req, MPI_STATUS_IGNORE);

}








void Wtime(int  Nf, clock_t begtt)
{

  if(rank==0){
time_t now;
time(&now);
  int hours1, minutes1, seconds1, day1, month1, year1;
struct tm *local = localtime(&now);
    hours1 = local->tm_hour;         
    minutes1 = local->tm_min;        
    seconds1 = local->tm_sec;        
    day1 = local->tm_mday;            
    month1 = local->tm_mon + 1;      
    year1 = local->tm_year + 1900;   

long double time_spent = 0.0;

clock_t endtt;

 int hours  =(int)time_spent/3600;
time_spent =time_spent-hours*3600;
int minutes=(int)(int)time_spent/60;
time_spent =time_spent-minutes*60;
int seconds=(int)time_spent;




    endtt = clock();
    time_spent = (long double)(endtt - begtt) / CLOCKS_PER_SEC;
    hours  =(int)time_spent/3600;
    time_spent =time_spent-hours*3600;
    minutes=(int)(int)time_spent/60;
    time_spent =time_spent-minutes*60;
    seconds=(int)time_spent;


    printf("%3d   %3d:%2d:%2d",Nf, hours,minutes,seconds);
    printf("      %02d/%02d/%d", day1, month1, (year1-2000));
    printf("      %02d:%02d:%02d ", hours1, minutes1, seconds1);
    printf(" \n");
   }
  }





  void Time_test( complex long double *Nun,complex long double *Nvn,complex long double *Nwn,complex long double *un,complex long double *vn,complex long double *wn,complex long double *p,complex long double * bm0 ,complex long double * bm1,complex long double * bm2,complex long double * Wb0,complex long double * Wb1,complex long double * Wb2){

clock_t begin = clock();


IC(Nun,Nvn,Nwn,un,vn,wn,p,bm0,bm1,bm2,Wb0,Wb1,Wb2);

for ( int nn = 0; nn <= 10; ++nn){

NSE_solver(Nun,Nvn,Nwn,un,vn,wn,p,bm0,bm1,bm2,Wb0,Wb1,Wb2);






}
clock_t end = clock();
clock_t begtt = clock();
    fft2mat(un);fft2mat(vn);fft2mat(wn);
    fft2mat(p);
    matrix_save(un,vn,wn,p,bm0,0);
    mat2fft(un);mat2fft(vn);mat2fft(wn);
    mat2fft(p);
    clock_t endtt = clock();



if(rank==0){

long double time_spent = 0.0;
time_spent += (N/10.0)*(long double)(end - begin) / CLOCKS_PER_SEC+(N)*(long double)(endtt - begtt) / CLOCKS_PER_SEC;

int hours  =(int)time_spent/3600;
time_spent =time_spent-hours*3600;
int minutes=(int)(int)time_spent/60;
time_spent =time_spent-minutes*60;
int seconds=(int)time_spent;

printf("The run time is:   %d:%d:%d \n", hours,minutes,seconds);
}


}




