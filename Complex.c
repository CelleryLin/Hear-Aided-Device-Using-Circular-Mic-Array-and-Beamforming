/*
gcc -Wall filename.c -lgsl -lgslcblas -lm
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>

#define ELEMENT 6
#define SAMPLERATE 44100
#define DEG2RAD(x) (x*M_PI/180)

double Radii=46.3e-3;
double LDA=340.f/1000.f;

void PrintMat(gsl_matrix *m){
    printf("Print Matrix: \n");
    for(int i=0; i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            printf(j==m->size2-1?"%6.2f\n":"%6.2f\t", gsl_matrix_get(m,i,j));
        }
    }
}

void SetMat(gsl_matrix *m, double n[]){
    //if( m->size1*m->size2 != (int)(sizeof(n)/sizeof(n[0]))){
    //    printf("Error: Set Matrix Mismatchd! gsl_len is %d and input array is %d.\n",
    //        m->size1*m->size2,
    //        (int)(sizeof(n)/sizeof(n[0])));
    //    exit(0);
    //}
    for(int i=0;i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            gsl_matrix_set (m, i, j, *n++);
        }
    }
}

void Beamforming(){
    double testsig[] = {1,2,3,4,5,
                        2,4,6,8,15,
                        3,6,10,12,40,
                        4,8,12,16,20,
                        5,10,14,20,25,
                        6,1,18,24,30};

    gsl_matrix *m=gsl_matrix_alloc(6, 5);
    gsl_matrix *mH=gsl_matrix_alloc(5, 6);
    gsl_matrix *Rxx=gsl_matrix_alloc(6, 6);

    SetMat(m,testsig);
    gsl_matrix_transpose_memcpy(mH, m);

    for (int i=0;i<m->size1;i++){
        for(int j=0;j<mH->size2;j++){
            gsl_vector_view tempa = gsl_matrix_row (m, i);
            gsl_vector_view tempb = gsl_matrix_column (mH, j);
            double corr = gsl_stats_correlation(tempa.vector.data, tempa.vector.stride,
                                                tempb.vector.data, tempb.vector.stride,
                                                tempa.vector.size);
            //printf("%f\n",corr);
            gsl_matrix_set (Rxx, i, j, corr);
        }
    }
    PrintMat(Rxx);
    
}


gsl_complex eulers_formula(double x){    //e^(xi) = cos(x) + i*sin(x)
    gsl_complex a;
    GSL_SET_COMPLEX(&a, cos(x), sin(x));
    return a;
}

void manifold(gsl_complex AA[], double doa){
    double phi = DEG2RAD(90);
    double temp;
    for(int k=0;k<ELEMENT;k++){
        temp=sin(phi)*cos(DEG2RAD(doa)-2*k*M_PI/ELEMENT);
        AA[k]=eulers_formula(-1*2*M_PI*Radii*temp/LDA);
        //printf("%f + %fi\n",AA[k]);
    }
}

int main(){
    gsl_complex AA[ELEMENT];
    manifold(AA, 0.);
    //for(int i=0;i<ELEMENT;i++)
    //    printf("%f+%fi\n",AA[i]);
    Beamforming();
    return 0;
}
