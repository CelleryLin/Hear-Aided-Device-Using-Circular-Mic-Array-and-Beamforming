/*
    gcc fftw_test.c -lfftw3 -lm
*/


#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <complex.h>
#include <math.h>

int main() {
    int M = 6;  // number of rows
    int N = 256;  // number of columns
    //double **input;  // input array
    //fftw_complex **output;  // output array (complex-valued)
    //fftw_plan *plan;  // plan for executing the FFT

    double input[M][N];
    fftw_complex output[M][N];
    fftw_plan plan[M];

    // Allocate memory for the input and output arrays
    //input = (double**) fftw_malloc(sizeof(double*) * M);
    //output = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*) * M);
    //plan = (fftw_plan*) fftw_malloc(sizeof(fftw_plan) * M);

    //for (int m = 0; m < M; m++) {
    //    input[m] = (double*) fftw_malloc(sizeof(double) * N);
    //    output[m] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
    //}

    // Initialize the input array

    double fs = 44100;  // sample rate (in Hz)
    double f = 440;  // frequency (in Hz)

    // Generate the sinusoidal signal
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            input[m][n] = sin(2*M_PI*f*n/fs);
        }
    }

    // Create a plan and executing the FFT
    for(int i=0;i<M;i++){
        plan[i] = fftw_plan_dft_r2c_1d(N, input[i], output[i], FFTW_ESTIMATE);
        fftw_execute(plan[i]);
    }


    // Print the output array
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N/2+1; n++) {
            printf("%d %d %f %f\n", m, n, output[m][n][0], output[m][n][1]);
        }
    }
}