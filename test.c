#include <stdio.h>
#include <math.h>
#include "portaudio.h"

int main(){
    int a[20]={1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,0,0};
    int bl[20]; int br[20];
    int *ap=a;
    int *blp=bl;
    int *brp=br;

    for(int i=0; i < 10; i++ ) {
        *blp++ = *ap++;
        *brp++ = *ap++;
        printf("%d, %d\n", bl[i],br[i]);
    }
    return 0;


}