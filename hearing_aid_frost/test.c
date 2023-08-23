/*
    gcc ./test.c -o ./a.out -lm
*/

# include <stdio.h>
# include <string.h>
# include <math.h>

#define FRAME_BLOCK_LEN 10
#define SAMPLING_RATE 48000
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1

# define M_PI 3.14159265358979323846
#define Radii 0.0463
#define SPEED 340

#define DOA 0 //rad
#define DELAY_STAGE_tau(n,doa) ((Radii/SPEED)*(cos(doa)-cos(doa-2*M_PI*n/USED_CH)))

float u[10] = { [0 ... 9-1] = 1. };
void matrixtest(float *a, int k, float n);

int main(){
    float a[3][1] = {0};
    // matrixtest(&a[0][0], 2, 3);
    // 
    // for(int i = 0; i < 3; i++){
    //     printf("%f\n", a[i][0]);
    // }
}

void matrixtest(float *a, int k, float n){
    a[k] = n;
}