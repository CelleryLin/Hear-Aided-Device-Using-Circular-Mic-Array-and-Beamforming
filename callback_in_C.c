/*
This program monitors audio in real-time.
Compile and link:

    gcc callback_in_C.c -lportaudio -lm

If use gsl:

    gcc callback_in_C.c -lportaudio -lgsl -lgslcblas -lm
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

#include "portaudio.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define FRAME_BLOCK_LEN 512
#define SAMPLING_RATE 44100
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1

PaStream *audioStream;

void SetMat(gsl_matrix_complex *m, gsl_complex n[]){
    for(int i=0;i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            gsl_matrix_complex_set (m, i, j, *n++);
        }
    }
}

gsl_complex double2complex(double re, double im){
    gsl_complex b;
    b.dat[0]=re;
    b.dat[1]=im;
    return b;
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


int audio_callback(const void *inputBuffer,
                   void *outputBuffer, 
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo* timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData){
    
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    static double phase = 0;
    float in_sort[INPUT_CHANNEL][FRAME_BLOCK_LEN];
    double fft_data[INPUT_CHANNEL][2*FRAME_BLOCK_LEN];
    double ifft_data[INPUT_CHANNEL][2*FRAME_BLOCK_LEN];
    gsl_matrix_complex *fft_data_mat=gsl_matrix_complex_alloc(INPUT_CHANNEL, FRAME_BLOCK_LEN);


    for (int i=0;i<INPUT_CHANNEL;i++){
        for (int j=0;j<FRAME_BLOCK_LEN;j++){
            in_sort[i][j]=in[INPUT_CHANNEL*j+i];
            REAL(fft_data[i],j) = in_sort[i][j];
            IMAG(fft_data[i],j) = 0.0;
        }
        gsl_fft_complex_radix2_forward(fft_data[i], 1, FRAME_BLOCK_LEN);

        for(int k=0;k<FRAME_BLOCK_LEN;k++){
            gsl_matrix_complex_set (fft_data_mat, i, k, double2complex(REAL(fft_data[i],k), IMAG(fft_data[i],k)));
        }
    }
    for (int i=0;i<INPUT_CHANNEL;i++){
        _gsl_vector_complex_view temp_view=gsl_matrix_complex_row(fft_data_mat,i);
        for (int j=0;j<FRAME_BLOCK_LEN;j++){
            gsl_complex temp_vew_obj=gsl_vector_complex_get((&temp_view.vector),j);
            REAL(ifft_data[i],j)=temp_vew_obj.dat[0];
            IMAG(ifft_data[i],j)=temp_vew_obj.dat[1];
        }

        gsl_fft_complex_radix2_inverse(ifft_data[i], 1, FRAME_BLOCK_LEN);
    }

    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        float temp=0;
        
        for(int j=0;j<INPUT_CHANNEL;j++){
            //temp+=*in++/6;
            temp+=in_sort[j][i]/6;
            //temp+=REAL(ifft_data[j],i)/6;

        }
        

        *out++ = temp;
        //*out++ = 0;
    }
    return paContinue;
}

void init_stuff(){
    int i,id;
    const PaDeviceInfo *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;

    Pa_Initialize(); // initialize portaudio

    id = Pa_GetDefaultOutputDevice();
    info = Pa_GetDeviceInfo(id);
    hostapi = Pa_GetHostApiInfo(info->hostApi);
    printf("output device [%s] %s\n",hostapi->name, info->name);
    outputParameters.device = id;
    outputParameters.channelCount = OUTPUT_CHANNEL;
    outputParameters.sampleFormat = paFloat32;
    outputParameters.suggestedLatency = 0;
    outputParameters.hostApiSpecificStreamInfo = NULL; /*no specific info*/

    sleep(0.5);
    id = Pa_GetDefaultInputDevice();
    info = Pa_GetDeviceInfo(id); /* get chosen device information struct */
    hostapi = Pa_GetHostApiInfo(info->hostApi); /* get host API struct */
    printf("input device [%s] %s\n",
    hostapi->name, info->name);
    inputParameters.device = id; /* chosen device id */
    inputParameters.channelCount = INPUT_CHANNEL; /* stereo input */
    inputParameters.sampleFormat = paFloat32; /* 32 bit float input */
    inputParameters.suggestedLatency = 0;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    Pa_OpenStream(
        &audioStream,
        &inputParameters,       // output parameters
        &outputParameters,      // input parameters
        SAMPLING_RATE,          // set sampling rate
        FRAME_BLOCK_LEN,        // set frames per buffer
        paClipOff,              // set no clip
        audio_callback,         // callback function
        NULL );                 // provide no data for the callback
    
    Pa_StartStream(audioStream); /* start the callback mechanism */
    printf("running. . . press space bar and enter to exit\n");
    printf("Latency: %f\n",info->defaultLowInputLatency);
}

void terminate_stuff(){
    Pa_StopStream( audioStream ); /* stop the callback mechanism */
    Pa_CloseStream( audioStream ); /* destroy the audio stream object */
    Pa_Terminate(); /* terminate portaudio */
}

int main(){
    init_stuff();
    while(getchar() != ' ') Pa_Sleep(100);
    terminate_stuff();
    return 0;
}



/*
def MVDR(input_buffer,doa):
    b=input_buffer[0:ELEMENT]
    CovMat = de.corr_matrix_estimate(b.T, imp="fast")
    doa_processing=doa
    if (LA.det(CovMat)!=0 and len(doa_processing)!=0):        
        invR=LA.inv(CovMat)
        temp=np.sin(phi)*np.cos(np.radians(doa_processing[0])-2*np.arange(0,ELEMENT)*np.pi/ELEMENT)
        AA=np.exp(-1j*2*np.pi*Radii*temp/LDA)
        ww=np.matmul(np.matmul(AA,invR),AA.conj().transpose())
        #print("ww : {}".format(ww))
        w=np.matmul(invR,AA.conj().transpose())*(1/ww)
        
        return w
        #output_buffer=np.matmul(w.conj().transpose(),b)
        #print(output_buffer.shape)
        #return np.real(output_buffer)
*/