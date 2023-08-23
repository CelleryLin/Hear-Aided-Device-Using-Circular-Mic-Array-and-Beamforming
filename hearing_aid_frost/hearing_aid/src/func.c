#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_linalg.h>
#include <fftw3.h>
#include <complex.h>
#include "params.h"


gsl_complex double2complex(float re, float im){
    gsl_complex b;
    b.dat[0]=re;
    b.dat[1]=im;
    return b;
}

float hanning(int n, int N){
    double w;
    w=0.5*(1-cos(2*M_PI*n/(N-1)));
    return (float)w;
}

float EMA(float a, float a0, float alpha){
    float new_val = a0 + alpha * (a - a0);
    return new_val;
}

float energy_vad(fftwf_complex a[],float threshold){
    float energy=0.;
    for(int i=0;i<FRAME_BLOCK_LEN/2+1;i++){
        energy+=a[i][0]*a[i][0]+a[i][1]*a[i][1];
    }
    if(energy>threshold){
        return 1.;
    }
    else{
        return 1e-9;    //prevent zero
    }
}

float SNR(fftwf_complex a[],double freq1, double freq2){
    float noise_avg=0.;
    float signal_avg=0.;
    float power_spectrum[FRAME_BLOCK_LEN/2+1];

    int start_n = (int)(freq1*FRAME_BLOCK_LEN/SAMPLING_RATE/2+1);  // freq to bins
    int end_n = (int)(freq2*FRAME_BLOCK_LEN/SAMPLING_RATE/2+1);

    for (int i=0; i<FRAME_BLOCK_LEN/2+1; i++) {
        power_spectrum[i] = a[i][0] * a[i][0] + a[i][1] * a[i][1];
    }
    
    
    for(int i=0;i<FRAME_BLOCK_LEN/2+1;i++){
        noise_avg+=power_spectrum[i];
    }
    noise_avg/=(FRAME_BLOCK_LEN/2+1);

    for(int i=start_n;i<=end_n;i++){
        signal_avg+=power_spectrum[i];
    }
    //printf("%f\n", signal_avg);
    signal_avg/=(end_n-start_n);

    //if((signal_avg/noise_avg)>3){
    //    printf("%f\n",(signal_avg/noise_avg));
    //}
    
    //printf("%f\t%f\n",signal_avg,noise_avg);
    //printf("%f\n",signal_avg/noise_avg);
    return (signal_avg/noise_avg);
    
    //return (signal_avg);
}

void PrintMat(gsl_matrix_complex *m, char info[]){
    printf("Print Matrix %s: \n",info);
    for(int i=0; i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            printf("%.6f+%.2fi\t", gsl_matrix_complex_get(m,i,j).dat[0],gsl_matrix_complex_get(m,i,j).dat[1]);
        }
        printf("\n");
    }
}

void autocorr(gsl_matrix_complex *Rxx, gsl_matrix_complex *m){
    gsl_matrix_complex *mH=gsl_matrix_complex_alloc(m->size2, m->size1);
    gsl_matrix_complex *Imat=gsl_matrix_complex_alloc(m->size1, m->size1);
    gsl_matrix_complex_set_identity(Imat);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,double2complex(1.,0.),m,Imat,double2complex(0.,0.),mH);
    //gsl_matrix_transpose_memcpy(mH, m);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
                    double2complex(1./m->size2,0.), m, mH,
                    double2complex(0.,0.), Rxx);

    gsl_matrix_complex_free(mH);
    gsl_matrix_complex_free(Imat);

}

void invMat(gsl_matrix_complex *matrix,gsl_matrix_complex *inv){
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);
    int s;

    gsl_linalg_complex_LU_decomp(matrix, p, &s);
    gsl_linalg_complex_LU_invert(matrix, p, inv);
    gsl_permutation_free(p);
}

gsl_complex eulers_formula(double x){    //e^(xi) = cos(x) + i*sin(x)
    gsl_complex a;
    GSL_SET_COMPLEX(&a, cos(x), sin(x));
    return a;
}
float Compress(float n){
    float a, b, tmp;
    a=Ka*n;
    b=Kb*n+1;
    return a/b;
}
void Compress_Mat(gsl_matrix_complex *matrix){
    gsl_complex a, b, tmp;
    
    for(int i=0; i<matrix->size1; i++){
        for(int j=0; j<matrix->size2; j++){
            if (gsl_matrix_complex_get(matrix, i, j).dat[0]<=Ka){
                // do nothing
            }
            else{
                gsl_matrix_complex_set(matrix, i, j, double2complex(1., 0.));
                //a = gsl_complex_mul(double2complex(Ka, 0.), gsl_matrix_complex_get(matrix, i, j));
                //b = gsl_complex_add(
                //gsl_complex_mul(
                //    double2complex(Kb, 0.),
                //    gsl_matrix_complex_get(matrix, i, j)),
                //double2complex(1., 0.));
                //tmp=gsl_complex_div(a, b);
                //gsl_matrix_complex_set(matrix, i, j, tmp);
            }
        }
    }
}