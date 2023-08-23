/*
This program monitors audio in real-time.
Compile and link:
    gcc ./src/main.c -Wall -Wextra -g -Iinclude -o ./a.out -lportaudio -lgsl -lgslcblas -lfftw3f -lm -Llib -Lsrc

*/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "params.h"
#include "matrix.h"
#include "portaudio.h"

float w[FILTER_LEN*USED_CH][1];
float F[FILTER_LEN*USED_CH][1];
float P[FILTER_LEN*USED_CH][FILTER_LEN*USED_CH];
float C[FILTER_LEN*USED_CH][FILTER_LEN];
float u[FILTER_LEN][1];


PaStream *audioStream;

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


	/* INPUT AND SORT */
    int ch_tag = 0;
    int zero_count = 0;
    for (int i = 0; i < INPUT_CHANNEL; i++)
    {
        float is_zero = 0.;
        for (int j = 0; j < FRAME_BLOCK_LEN; j++)
        {
            in_unsort[i][j] = in[INPUT_CHANNEL * j + i];
            is_zero += (float)fabs((double)in[INPUT_CHANNEL * j + i]);
        }
        if (is_zero < 1e-5 && zero_count != 2)
        {
            zero_count++;
            continue;
        }
        if (zero_count >= 2)
        {
            for (int j = 0; j < FRAME_BLOCK_LEN; j++)
            {
                in_sort[ch_tag][j] = in_unsort[i][j];
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

    /* STEERING */
    float steering_output_T[USED_CH*FILTER_LEN][FRAME_BLOCK_LEN];
    for (int i=0;i<USED_CH;i++){
        for (int j=0;j<FILTER_LEN;j++){
            int sample_shift = (int)(DELAY_STAGE_tau(i,(double)(DOA)*M_PI/180)*SAMPLING_RATE);
            signal_shift(in_sort[i], steering_output_T[i+j*USED_CH], sample_shift+j); 
        }
    }

	/* OUTPUT */
    float original[FRAME_BLOCK_LEN];
    float *p_original = original;

    float filtered[1][FRAME_BLOCK_LEN];
    float *p_filtered = filtered[0];

    // check if w is too big
    for(int i=0;i<FILTER_LEN*USED_CH;i++){
        if(fabs(w[i][0])>10){
            // set w to F
            for(int i=0;i<FILTER_LEN*USED_CH;i++){
                w[i][0]=F[i][0];
            }
            break;
        }
    }

    FIR_filter(&steering_output_T[0][0], &filtered[0][0], &w[0][0]);

    for (int i=0;i<36;i++){
       printf("%f ",w[i][0]);
    }
    printf("\n");

    for (int i = 0; i < FRAME_BLOCK_LEN; i++){
        // float tmp=0;
        // for(int j=0;j<USED_CH;j++){
        //     tmp+=(steering_output_T[j][i]/6);
        // }
        // original[i]=tmp;
        *out++ = *p_filtered++;
        //*out++ = *p_steering_output++;
    }

    /* FROST BEAMFORMING */

    float steering_output[FRAME_BLOCK_LEN][USED_CH*FILTER_LEN];
    transpose_matrix(&steering_output_T[0][0], &steering_output[0][0], USED_CH*FILTER_LEN, FRAME_BLOCK_LEN);
    float yX[1][USED_CH*FILTER_LEN];
    matmul(&filtered[0][0], &steering_output[0][0], &yX[0][0], 1, FRAME_BLOCK_LEN, USED_CH*FILTER_LEN);

    float yX_T[USED_CH*FILTER_LEN][1];
    transpose_matrix(&yX[0][0], &yX_T[0][0], 1, USED_CH*FILTER_LEN);

    for(int i=0;i<USED_CH*FILTER_LEN;i++){
        yX_T[0][i]*=LEARNING_RATE;
    }

    float WmuyX[USED_CH*FILTER_LEN][1];
    matsub(&w[0][0], &yX_T[0][0], &WmuyX[0][0],  1,USED_CH*FILTER_LEN);

    float PWmuyX[USED_CH*FILTER_LEN][1];
    matmul(&P[0][0], &WmuyX[0][0], &PWmuyX[0][0], USED_CH*FILTER_LEN, USED_CH*FILTER_LEN, 1);


    float PWmuyXF[USED_CH*FILTER_LEN][1];
    matadd(&PWmuyX[0][0], &F[0][0], &PWmuyXF[0][0], USED_CH*FILTER_LEN, 1);

    memcpy(&w[0][0], &PWmuyXF[0][0], sizeof(float)*USED_CH*FILTER_LEN);


    return paContinue;
}

void init_stuff()
{
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


    /* INIT BEAMFORM PATAMS*/
    // u = ones
    // C = [c1, c2, c3, ... cFILTER_LEN], cj = [0,0,0,...,0,1,1,1,...1,0,0,0,...0]' contains USED_CH 1s at jth group
    // w = w_i = C*inv(C'*C)*u = F
    // P = I - C*inv(C'*C)*C'

    //u
    set_u(&u[0][0]);

    // C
    float C_T[FILTER_LEN][FILTER_LEN*USED_CH] = {{0.}};
    float ones[USED_CH] = {[0 ... USED_CH-1] = 1.};
    for (int i=0;i<FILTER_LEN;i++){
        memcpy(&C_T[i][i*USED_CH],ones,sizeof(ones[0])*USED_CH);
    }

    transpose_matrix(&C_T[0][0],&C[0][0],FILTER_LEN,FILTER_LEN*USED_CH);

    // // print C
    // for(int i=0;i<FILTER_LEN*USED_CH;i++){
    //     for(int j=0;j<FILTER_LEN;j++){
    //         printf("%f ",C[i][j]);
    //     }
    //     printf("\n");
    // }

    // F
    float I[FILTER_LEN*USED_CH][FILTER_LEN*USED_CH] = {{0.}};
    for (int i=0;i<FILTER_LEN*USED_CH;i++){
        I[i][i]=1.;
    }

    float CTC[FILTER_LEN][FILTER_LEN];
    matmul(&C_T[0][0],&C[0][0],&CTC[0][0],FILTER_LEN,FILTER_LEN*USED_CH,FILTER_LEN);

    float invCTC[FILTER_LEN][FILTER_LEN];
    calculateInverse(&CTC[0][0],&invCTC[0][0],FILTER_LEN);

    float CinvCTC[FILTER_LEN*USED_CH][FILTER_LEN];
    matmul(&C[0][0],&invCTC[0][0],&CinvCTC[0][0],FILTER_LEN*USED_CH,FILTER_LEN,FILTER_LEN);

    matmul(&CinvCTC[0][0],&u[0][0],&F[0][0],FILTER_LEN*USED_CH,FILTER_LEN,1);

    for(int i=0;i<FILTER_LEN*USED_CH;i++){
        w[i][0]=F[i][0];
    }

    // P
    float CinvCTCCT[FILTER_LEN*USED_CH][FILTER_LEN*USED_CH];
    
    matmul(&CinvCTC[0][0],&C_T[0][0],&CinvCTCCT[0][0],FILTER_LEN*USED_CH,FILTER_LEN,FILTER_LEN*USED_CH);
    matsub(&I[0][0],&CinvCTCCT[0][0],&P[0][0],FILTER_LEN*USED_CH,FILTER_LEN*USED_CH);


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

void signal_shift(float *input, float *output, int tau){
    // tau is the shift amount (sampling points)
    // output(k) = input(k+tau)

    // float tmp[FRAME_BLOCK_LEN] = {0.};
    if (tau>=0){
        memcpy(output, &input[tau], sizeof(output)*(FRAME_BLOCK_LEN - tau));

        //wrap
        //memcpy(&output[FRAME_BLOCK_LEN-tau], input, sizeof(output)*(tau));
        
        //zero padding
        //memcpy(&output[FRAME_BLOCK_LEN-tau], tmp, sizeof(output)*(tau));

        //padding from original
        memcpy(&output[FRAME_BLOCK_LEN-tau], &input[FRAME_BLOCK_LEN-tau], sizeof(output)*(tau));

        //morphing (interpolation)
        if(FRAME_BLOCK_LEN-tau+1 <= FRAME_BLOCK_LEN){
            float diff = output[FRAME_BLOCK_LEN-tau+1] - output[FRAME_BLOCK_LEN-tau-2];
            output[FRAME_BLOCK_LEN-tau-1] = output[FRAME_BLOCK_LEN-tau-2] + diff*0.333;
            output[FRAME_BLOCK_LEN-tau] = output[FRAME_BLOCK_LEN-tau-2] + diff*0.666;
        }
    }
    else{
        tau = -tau;
        memcpy(&output[tau], input, sizeof(output)*(FRAME_BLOCK_LEN - tau));

        //wrap
        //memcpy(output, &input[FRAME_BLOCK_LEN-tau], sizeof(output)*(tau));
        
        //zero padding
        //memcpy(output, &tmp[FRAME_BLOCK_LEN-tau], sizeof(output)*(tau));

        //padding from original
        memcpy(output, &input, sizeof(output)*(tau));

        //morphing (interpolation)
        if(tau-2 >= 0){
            float diff = output[tau-2] - output[tau+1];
            output[tau-1] = output[tau+1] + diff*0.333;
            output[tau] = output[tau+1] + diff*0.666;
        }
    }
}

void FIR_filter(float *input, float *output, float *w){
    /*
    input = [[s1mic1],[s1mic2],...,[s1micN],[s2mic1],[s2mic2],...,[s2micN],...,[sLmicN]]
    which N is the number of microphones and L is the number of Filter taps
    input size is [N*L][FRAME_BLOCK_LEN] 2D matrix

    w = [w1,w2,...,wN,wN+1,wN+2...,wN*L]'
    w size is [N*L] column vector

    output is the filter output, size is [FRAME_BLOCK_LEN] row vector
    */
    
    for (int i=0;i<FRAME_BLOCK_LEN;i++){
        float tmp = 0;
        for (int j=0;j<FILTER_LEN*USED_CH;j++){
            tmp += input[i+j*FRAME_BLOCK_LEN]*w[j];
        }
        output[i] = tmp;
    }
}

void set_u(float *fir){
    fir[0] = 1.;
    fir[1] = 1.;
    fir[2] = 1.;
    fir[3] = 1.;
    fir[4] = 1.;
    fir[5] = 1.;
}

int main(){
    init_stuff();
    while (getchar() != ' ')
        Pa_Sleep(100);
    terminate_stuff();
    return 0;
}