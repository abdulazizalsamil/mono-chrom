#include"include/main.h"

double kappa(int j){
double w;

if (j<Ypt/2+1){w= -j/Ly;}
else          {w= (Ypt-j)/Ly;}


  return w;
}

/*(1)
 the function (complex2FFT) is taking a vector U as a 
 complex (actually as reall) and return the FFT
 conversion.

*/
void complex2FFT( complex long double *U)
{


/*
complex2FFT is used to convert a complex vector U to FFT:
First, defining an (in), and (out) fftwl_complex vectors then, assiging 
the values from U to (in). The fft_plan is used to convert the (in) to 
fft. the values of (out) are assigned back to U and the normalizing_factor
is used to normalize the data. 
 Note: the letter l in fftwl_complex revere to long double type of data. 
*/

      fftwl_complex *in;
      fftwl_complex *out;

  in  = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  out = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  
  for (int i = 0; i < Ypt; ++i){
      in[i ][0]= creal(U[i]);
      in[i ][1]= cimag(U[i]);
      }
fft(in, out);

  for (int i = 0; i < Ypt; ++i){U[i] = out[i][0] + out[i][1]*I;}

  fftwl_free(in);
  fftwl_free(out);


}


/*(2)
 the function (FFT2complex) is taking a FFT conversion of U 
 as a complex  and return a complex vector (actually as reall)
*/

void FFT2complex( complex long double *U)
{

/*
FFT2complex is used to convert FFT to  complex vector U:First, defining 
an (in), and (out) fftwl_complex vectors then, assiging the values from 
U to (in). The fft_plan is used to convert the (in) to complex vector. 
the values of (out) are assigned back to U and the normalizing_factor
is used to normalize the data. 
 Note: the letter l in fftwl_complex revere to long double type of data. 

FFT_test: The returned value (in) should contain zeros in the imaginary 
part.
*/
    
      fftwl_complex *in;
      fftwl_complex *out;

  in  = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  out = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  
  for (int i = 0; i < Ypt; ++i){
      out[i ][0]= creal(U[i]);
      out[i ][1]= cimag(U[i]);
      }

 
 
ifft(in,out);

for (int i = 0; i < Ypt; ++i){ U[i] =in[i ][0] + in[i ][1]*I;}


  fftwl_free(in);
  fftwl_free(out);
}



/*(3)
The functions mat2fft and fft2mat take matrix and covert to FFT and back
to complex matrix
*/


void mat2fft( complex long double *un){

/*
  fftwl_complex *in;
  fftwl_complex *out;
  in  = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  out = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
  fftwl_plan plan;


for (int i = ib-1; i < ie+1; ++i){
  for (int k = 0; k < Zpt; ++k){
  
    for (int j = 0; j < Ypt; ++j){
       in[j ][0]= creal(Un(i,j,k));
       in[j ][1]= cimag(Un(i,j,k));
       }

   plan =fftwl_plan_dft_1d(Ypt, in, out, FFTW_FORWARD , FFTW_ESTIMATE);
   fftwl_execute(plan);

    for (int j = 0; j < Ypt; ++j){Un(i,j,k) = (out[j][0] + out[j][1]*I)/(Ypt);}
    
     }}

  fftwl_free(in);
  fftwl_free(out);
  fftwl_destroy_plan(plan);
  fftwl_cleanup();
*/

}



void fft2mat( complex long double *un){

/*
 fftwl_complex *in;
 fftwl_complex *out;
 in  = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
 out = (fftwl_complex *) fftwl_malloc(sizeof(fftwl_complex) * Ypt);
 fftwl_plan plan;

for (int i = ib-1; i < ie+1; ++i){
  for (int k = 0; k < Zpt; ++k){


  
   for (int j = 0; j < Ypt; ++j){ 
     out[j ][0]= creal(Un(i,j,k));
     out[j ][1]= cimag(Un(i,j,k));
    }

 
    plan =fftwl_plan_dft_1d(Ypt, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftwl_execute(plan);


    for (int j = 0; j < Ypt; ++j){ 
     Un(i,j,k)=in[j ][0] + in[j ][1]*I;
    }


    }}


  fftwl_free(in);
  fftwl_free(out);
  fftwl_destroy_plan(plan);
  fftwl_cleanup();
*/
}


/*(4)
 the function (FFT_test) is used to test if the FFTW
 working by sending an sine valuses of Y to FFT and 
 inveversing it back to the original values.  

*/

void FFT_test( )
{
	/*
FFT_test is used to test the FFT: Two complex varible (U,V) assigned 
the same values and convert the U to complex and back to complex. Then,
the values of U and U are compared. The accuracy is of order 10^(-16).

	*/

complex long double *U = malloc(Ypt * sizeof(*U)); 
complex long double *V = malloc(Ypt * sizeof(*V)); 

  for (int i = 0; i < Ypt; ++i){
      U[i] = sin( (double) (2.0*PI* (double)i *Ly* dy)) ;
      V[i] = sin( (double) (2.0*PI* (double)i *Ly* dy)) ;


    }



complex2FFT(U);


FFT2complex(U);


for (int i = 0; i < Ypt; ++i){
   if(( fabs( creal(U[i]-V[i]) ) >pow(10,-16) ||  fabs( cimag(U[i]-V[i]))>pow(10,-16))){
  printf("FFTW_ERRO: The FFTW3 library is not is not working correctly.\n");
  printf("           Check functions (complex2FFT) and (FFT2complex) in fftw_f.c\n");
   exit(1); }}






free(U);
free(V);

}




 void fft(fftwl_complex *in, fftwl_complex *out)
{

  fftwl_plan plan =fftwl_plan_dft_1d(Ypt, in, out, FFTW_FORWARD , FFTW_ESTIMATE);

  fftwl_execute(plan);
  fftwl_destroy_plan(plan);
  fftwl_cleanup();
    for (int i = 0; i < Ypt; ++i) {
    out[i][0] /= (Ypt);
    out[i][1] /= (Ypt);
  }
}


void ifft(fftwl_complex *in, fftwl_complex *out)
{

  fftwl_plan plan =fftwl_plan_dft_1d(Ypt, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwl_execute(plan);

  fftwl_destroy_plan(plan);
  fftwl_cleanup();
}





