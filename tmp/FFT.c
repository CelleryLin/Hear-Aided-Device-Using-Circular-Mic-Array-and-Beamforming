/*
gcc -Wall filename.c -lgsl -lgslcblas -lm
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define BUFFER_SIZE 128
int main() {
    int i; double data[2*BUFFER_SIZE];

    for (i = 0; i < BUFFER_SIZE; i++){
    REAL(data,i) = 0.0; IMAG(data,i) = 0.0;
    }

    REAL(data,0) = 1.0;

    for (i = 1; i <= 10; i++){
        REAL(data,i) = REAL(data,BUFFER_SIZE-i) = 1.0;
    }
    
    for (i = 0; i < BUFFER_SIZE; i++){
        printf ("%d %e %e\n", i, REAL(data,i), IMAG(data,i));
    }

    printf ("\n\n");

    gsl_fft_complex_radix2_forward (data, 1, BUFFER_SIZE);

    for (i = 0; i < BUFFER_SIZE; i++){
        printf ("%d %e %e\n", i,
            REAL(data,i)/sqrt(BUFFER_SIZE),
            IMAG(data,i)/sqrt(BUFFER_SIZE));
    }

  return 0;
}