/*
gcc callback_in_C.c -lportaudio -lm
*/

#include <stdio.h>
#include <math.h>
#include "portaudio.h"

#define FRAME_BLOCK_LEN 256
#define SAMPLING_RATE 44100
#define TWO_PI (3.14159265f * 2.0f)

PaStream *audioStream;
double si = 0;

int main(){
    int id;
    const PaDeviceInfo *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    Pa_Initialize();

    for (int i=0;i < Pa_GetDeviceCount(); i++) {
        info = Pa_GetDeviceInfo(i); /* get information from current device */
        hostapi = Pa_GetHostApiInfo(info->hostApi); /*get info from curr. host API */
        if (info->maxOutputChannels > 0)/* if curr device supports output */
            printf("%d: [%s] %s (output)\n",i, hostapi->name, info->name );
    }
    id = Pa_GetDefaultDevice();
    printf("%d\n",id);
}