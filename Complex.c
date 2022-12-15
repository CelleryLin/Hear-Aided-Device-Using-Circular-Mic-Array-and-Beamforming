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
#include <gsl/gsl_linalg.h>

#define ELEMENT 6
#define SAMPLERATE 44100
#define DEG2RAD(x) (x*M_PI/180)

double Radii=46.3e-3;
double LDA=340.f/1000.f;

gsl_complex double2complex(double a){
    gsl_complex b;
    b.dat[0]=a;
    b.dat[1]=0.;
    return b;
}

void invMat(gsl_matrix_complex *matrix,gsl_matrix_complex *inv){
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);
    //gsl_matrix *mattemp=gsl_matrix_alloc(matrix->size1, matrix->size2);;
    //gsl_matrix_memcpy(mattemp,matrix);
    int s;

    gsl_linalg_complex_LU_decomp(matrix, p, &s);
    gsl_linalg_complex_LU_invert(matrix, p, inv);
    gsl_permutation_free(p);
}

void autocorr(gsl_matrix_complex *Rxx, gsl_matrix_complex *m){
    gsl_matrix_complex *mH=gsl_matrix_complex_alloc(m->size2, m->size1);
    gsl_matrix_complex *Imat=gsl_matrix_complex_alloc(m->size1, m->size1);
    gsl_matrix_complex_set_identity(Imat);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,double2complex(1.),m,Imat,double2complex(0.),mH);
    //gsl_matrix_transpose_memcpy(mH, m);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
                    double2complex(1./m->size2), m, mH,
                    double2complex(0.), Rxx);

    //for (int i=0;i<m->size1;i++){
    //    for(int j=0;j<mH->size2;j++){
    //        gsl_vector_view tempa = gsl_matrix_row (m, i);
    //        gsl_vector_view tempb = gsl_matrix_column (mH, j);
    //        double corr = gsl_stats_correlation(tempa.vector.data, tempa.vector.stride,
    //                                            tempb.vector.data, tempb.vector.stride,
    //                                            tempa.vector.size);
    //        //printf("%f\n",corr);
    //        gsl_matrix_set (Rxx, i, j, corr);
    //    }
    //}
}


void PrintMat(gsl_matrix_complex *m, char info[]){
    printf("Print Matrix %s: \n",info);
    for(int i=0; i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            printf("%6.2e+%.2fi\t", gsl_matrix_complex_get(m,i,j).dat[0],gsl_matrix_complex_get(m,i,j).dat[1]);
        }
        printf("\n");
    }
}

void SetMat(gsl_matrix_complex *m, gsl_complex n[]){
    for(int i=0;i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            gsl_matrix_complex_set (m, i, j, *n++);
        }
    }
}

void Beamforming(){
    double testsig[] = {1,2,3,4,5,
                        1.1,2,3.5,4,6,
                        2,3,4,5,1.1,
                        2,3.5,1,4,6,
                        3,4,5.5,6,1,
                        3.1,4,4,5,6};


    gsl_complex testsig_comp[30];
    for(int i=0;i<30;i++){
        testsig_comp[i]=double2complex(testsig[i]);
    }

    gsl_matrix_complex *m=gsl_matrix_complex_alloc(6, 5);
    gsl_matrix_complex *Rxx=gsl_matrix_complex_alloc(m->size1, m->size1);
    gsl_matrix_complex *invR=gsl_matrix_complex_alloc(m->size1, m->size1);
    SetMat(m,testsig_comp);
    // SetMat(m,testsig);
    autocorr(Rxx,m);
    //gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasNoTrans,CblasUnit,1.,Rxx,invR);
    PrintMat(Rxx,"Rxx");

    invMat(Rxx,invR);   // Note: this will change the val of Rxx!

    PrintMat(invR,"invR");
    
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
