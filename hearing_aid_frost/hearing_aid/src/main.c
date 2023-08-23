/*
This program monitors audio in real-time.
Compile and link:

    gcc callback_in_C.c -lportaudio -lm

If use gsl:

    gcc callback_in_C.c -lportaudio -lgsl -lgslcblas -lm

If use fftwf:
    gcc callback_in_C_wide.c -lportaudio -lgsl -lgslcblas -lfftw3f -lm
*/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
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
#include "portaudio.h"

#include "main.h"
#include "only_bf.h"
#include "only_bf_terminate.h"
#include "rt_nonfinite.h"


#define REAL(z, i) ((z)[2 * (i)])
#define IMAG(z, i) ((z)[2 * (i) + 1])
#define DEG2RAD(x) (x * M_PI / 180)

float snr = 1.;
float snr0 = 0.;
float snr_Intergral_val = 0.;
float snr_Previous_error = 0.;
float hanning_window[FRAME_BLOCK_LEN];
float doa = 1.;
float doa0 = 0.;

PaStream *audioStream;

/************************************************************************************/
/* Function Declarations */
static void argInit_1200x6_real_T(double result[7200]);

static double argInit_real_T(void);

/* Function Definitions */
static void argInit_1200x6_real_T(double result[7200])
{

    int i;
    for (i = 0; i < 7200; i++)
    {
        result[i] = argInit_real_T();
    }
}

static double argInit_real_T(void)
{
    return 0.0;
}

/************************************************************************************/


int audio_callback(const void *inputBuffer,
                   void *outputBuffer,
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData)
{

    float *in = (float *)inputBuffer, *out = (float *)outputBuffer;
    float in_unsort[INPUT_CHANNEL][FRAME_BLOCK_LEN];
    float in_sort[USED_CH][FRAME_BLOCK_LEN];

    /* SORT INPUT */
    int ch_tag = 0;
    int zero_count = 0;
    for (int i = 0; i < INPUT_CHANNEL; i++)
    {
        float is_zero = 0.;
        for (int j = 0; j < FRAME_BLOCK_LEN; j++)
        {
            in_unsort[i][j] = in[INPUT_CHANNEL * j + i];
            is_zero += (float)fabs((double)in[INPUT_CHANNEL * j + i]);
            // is_zero+=in[INPUT_CHANNEL*j+i];
            // printf("%f\t",in[INPUT_CHANNEL*j+i]);
        }
        // printf("\n\n");
        if (is_zero < 1e-5 && zero_count != 2)
        {
            // printf("WHAT THE HECK DID YOU DOl?");
            zero_count++;
            continue;
        }
        // printf("\n---DONE!---\n");
        if (zero_count >= 2)
        {
            for (int j = 0; j < FRAME_BLOCK_LEN; j++)
            {
                in_sort[ch_tag][j] = in_unsort[i][j];
                // in_sort_window[ch_tag][j] = in_unsort[i][j] * hanning_window[j];
            }
            // printf("\n\n");
            ch_tag++;
        }
    }
    for (int i = ch_tag; i < USED_CH; i++)
    {
        for (int j = 0; j < FRAME_BLOCK_LEN; j++)
        {
            in_sort[i][j] = in_unsort[i - ch_tag][j];
            // in_sort_window[i][j] = in_unsort[i - ch_tag][j] * hanning_window[j];
        }
        // printf("\n\n");
    }
    /* BEAMFORMING */
    double dv[7200];
    double output[1200];
    float f_output[1200];
    float *p_output = f_output;

    for (int i=0;i<USED_CH;i++){
        for (int j=0;j<FRAME_BLOCK_LEN;j++){
            dv[i*FRAME_BLOCK_LEN+j] = in_sort[i][j];
        }
    }

    argInit_1200x6_real_T(dv);
    only_bf(dv, output);

    for (int i=0;i<1200;i++){
        f_output[i] = (float)(output[i]);
    }

    /* OUTPUT */
    float original[FRAME_BLOCK_LEN];
    float *p_original = original;
    for (int i = 0; i < FRAME_BLOCK_LEN; i++){
        // float tmp=0;
        // for(int j=0;j<USED_CH;j++){
        //     tmp+=(in_sort[j][i]/6);
        //     //tmp+=output[j][i];
        // }
        // original[i]=tmp;
        *out++ = *p_output++;
    }
    return paContinue;
}

void init_stuff(){

    int id;
    const PaDeviceInfo *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    Pa_Initialize(); // initialize portaudio

    const PaDeviceInfo *deviceInfo;
    int numDevices;
    numDevices = Pa_GetDeviceCount();
    for (int i = 0; i < numDevices; i++)
    {
        deviceInfo = Pa_GetDeviceInfo(i);
        printf("%s\n", deviceInfo->name);
    }
    id = Pa_GetDefaultOutputDevice();
    info = Pa_GetDeviceInfo(id);
    hostapi = Pa_GetHostApiInfo(info->hostApi);
    printf("output device [%s] %s\n", hostapi->name, info->name);
    outputParameters.device = id;
    outputParameters.channelCount = OUTPUT_CHANNEL;
    outputParameters.sampleFormat = paFloat32;
    outputParameters.suggestedLatency = 0;
    outputParameters.hostApiSpecificStreamInfo = NULL; /*no specific info*/

    sleep(0.5);
    id = Pa_GetDefaultInputDevice();
    info = Pa_GetDeviceInfo(id);                /* get chosen device information struct */
    hostapi = Pa_GetHostApiInfo(info->hostApi); /* get host API struct */
    printf("input device [%s] %s\n",
           hostapi->name, info->name);
    inputParameters.device = id;                  /* chosen device id */
    inputParameters.channelCount = INPUT_CHANNEL; /* stereo input */
    inputParameters.sampleFormat = paFloat32;     /* 32 bit float input */
    inputParameters.suggestedLatency = 0;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    for (int i=0;i<FRAME_BLOCK_LEN;i++){
        hanning_window[i]=hanning(i,FRAME_BLOCK_LEN);
    }

    Pa_OpenStream(
        &audioStream,
        &inputParameters,  // output parameters
        &outputParameters, // input parameters
        SAMPLING_RATE,     // set sampling rate
        FRAME_BLOCK_LEN,   // set frames per buffer
        paClipOff,         // set no clip
        audio_callback,    // callback function
        NULL);             // provide no data for the callback

    Pa_StartStream(audioStream); /* start the callback mechanism */
    printf("running. . . press space bar and enter to exit\n");
    printf("Latency: %f\n", info->defaultLowInputLatency);
}

void terminate_stuff(){

    Pa_StopStream(audioStream);  /* stop the callback mechanism */
    Pa_CloseStream(audioStream); /* destroy the audio stream object */
    Pa_Terminate();              /* terminate portaudio */
}

int main(int argc, char **argv){

    init_stuff();
    while (getchar() != ' ')
        Pa_Sleep(100);
    terminate_stuff();
    return 0;
}
