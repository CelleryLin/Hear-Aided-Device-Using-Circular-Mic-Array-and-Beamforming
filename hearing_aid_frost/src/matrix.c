#include <stdio.h>

#include "matrix.h"
#include "params.h"


void transpose_matrix(float *input, float *output, int m, int n){
    // printf("m: %d, n: %d\n",m,n);
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            output[j*m+i] = input[i*n+j];
        }
    }
}

void matmul(float *inputA, float *inputB, float *output, int k, int m, int n){
    // inputA [k][m]
    // inputB [m][n]
    // output = A*B [k][n]
    for (int i=0;i<k;i++){
        for (int j=0;j<n;j++){
            float tmp=0;
            for (int l=0;l<m;l++){
                tmp+=inputA[i*m+l]*inputB[l*n+j];
            }
            output[i*n+j]=tmp;
        }
    }
}

void matsub(float *inputA, float *inputB, float *output, int m, int n){
    // inputA [m][n]
    // inputB [m][n]
    // output = A-B [m][n]
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            output[i*n+j]=inputA[i*n+j]-inputB[i*n+j];
        }
    }
}

void matadd(float *inputA, float *inputB, float *output, int m, int n){
    // inputA [m][n]
    // inputB [m][n]
    // output = A+B [m][n]
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            output[i*n+j]=inputA[i*n+j]+inputB[i*n+j];
        }
    }
}


float determinant(float *input, int k){
    float s = 1, det = 0, b[k-1][k-1];
    int m, n;
    if (k == 1){
        return (input[0]);
    }
    else{
        det = 0;
        for (int c = 0; c < k; c++){
            m = 0;
            n = 0;
            for (int i = 0; i < k; i++){
                for (int j = 0; j < k; j++){
                    // b[i][j] = 0;
                    if (i != 0 && j != c){
                        b[m][n] = input[i*k + j];
                        if (n < (k - 2))
                            n++;
                        else{
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (input[c] * determinant(&b[0][0], k - 1));
            s = -1 * s;
        }
    }
    return (det);
}

void calculateInverse(float *input, float *output, int k) {

    float det = determinant(input, k);
    if (det <= 1e-7 && det >= -1e-7) {
        printf("\nInverse of Entered Matrix is not possible\n");
        return;
    }

    float augmentedMatrix[k][2 * k];

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; ++j) {
            augmentedMatrix[i][j] = input[i*k + j];
            augmentedMatrix[i][j + k] = (i == j) ? 1.0f : 0.0f;
        }
    }

    for (int i = 0; i < k; i++) {
        float diagElement = augmentedMatrix[i][i];
        for (int j = 0; j < 2 * k; j++) {
            augmentedMatrix[i][j] /= diagElement;
        }

        for (int j = 0; j < k; j++) {
            if (j != i) {
                float factor = augmentedMatrix[j][i];
                for (int l = 0; l < 2 * k; l++) {
                    augmentedMatrix[j][l] -= factor * augmentedMatrix[i][l];
                }
            }
        }
    }

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            output[i*k + j] = augmentedMatrix[i][j + k];
        }
    }
}