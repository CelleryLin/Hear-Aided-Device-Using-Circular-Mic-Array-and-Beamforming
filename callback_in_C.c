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
#include <fftw3.h>

#include "portaudio.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define DEG2RAD(x) (x*M_PI/180)

#define FRAME_BLOCK_LEN 1024
#define SAMPLING_RATE 44100
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1

#define LDA 343./1000.
#define Radii 0.0463
#define Ksnr 0.05
#define GAIN 2

PaStream *audioStream;
double snr=1.;
    
gsl_complex double2complex(double re, double im){
    gsl_complex b;
    b.dat[0]=re;
    b.dat[1]=im;
    return b;
}

float SNR(double a[],double m_freq, double bandwidth){
    double noise_avg=0.;
    double signal_avg=0.;
    double power_spectrum[FRAME_BLOCK_LEN];
    int m_n = (int)(m_freq*FRAME_BLOCK_LEN/SAMPLING_RATE)+(FRAME_BLOCK_LEN/2)+1;
    int bandwidth_n = (int)(bandwidth*FRAME_BLOCK_LEN/SAMPLING_RATE);

    for (int i=0; i<FRAME_BLOCK_LEN; i++) {
        power_spectrum[i] = REAL(a,i) * REAL(a,i) + IMAG(a,i) * IMAG(a,i);
    }
    
    
    for(int i=(FRAME_BLOCK_LEN/2)+1;i<FRAME_BLOCK_LEN;i++){
        noise_avg+=fabs(REAL(power_spectrum,i));
    }
    noise_avg/=(FRAME_BLOCK_LEN/2);

    for(int i=m_n-bandwidth_n;i<=m_n+bandwidth_n;i++){
        signal_avg+=fabs(REAL(power_spectrum,i));
    }
    signal_avg/=bandwidth_n;

    //if((signal_avg/noise_avg)>3){
    //    printf("%f\n",(signal_avg/noise_avg));
    //}
    
    printf("%f\t%f\n",signal_avg,noise_avg);
    printf("%f\n",power_spectrum[100]);
    //return (signal_avg/noise_avg);
    //return (signal_avg);
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

void autocorr(gsl_matrix_complex *Rxx, gsl_matrix_complex *m){
    gsl_matrix_complex *mH=gsl_matrix_complex_alloc(m->size2, m->size1);
    gsl_matrix_complex *Imat=gsl_matrix_complex_alloc(m->size1, m->size1);
    gsl_matrix_complex_set_identity(Imat);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,double2complex(1.,0.),m,Imat,double2complex(0.,0.),mH);
    //gsl_matrix_transpose_memcpy(mH, m);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
                    double2complex(1./m->size2,0.), m, mH,
                    double2complex(0.,0.), Rxx);
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

gsl_complex eulers_formula(double x){    //e^(xi) = cos(x) + i*sin(x)
    gsl_complex a;
    GSL_SET_COMPLEX(&a, cos(x), sin(x));
    return a;
}


int audio_callback(const void *inputBuffer,
                   void *outputBuffer, 
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo* timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData){
    
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    double y_arr_f_to_t_domain[2*FRAME_BLOCK_LEN];
    double *y_arr_p = y_arr_f_to_t_domain;
    static double phase = 0;
    float in_unsort[INPUT_CHANNEL][FRAME_BLOCK_LEN];
    float in_sort[USED_CH][FRAME_BLOCK_LEN];
    double fft_data[USED_CH][2*FRAME_BLOCK_LEN];
    double ifft_data[USED_CH][2*FRAME_BLOCK_LEN];
    int printmat=0;
    gsl_matrix_complex *fft_data_mat=gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    gsl_matrix_complex *in_sort_mat=gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    plan = fftw_plan_dft_r2c_1d(FRAME_BLOCK_LEN, input, output, FFTW_ESTIMATE);
    //printf("\n\n-----------start print------------\n");
    int ch_tag=0;
    int zero_count=0;
    for (int i=0;i<INPUT_CHANNEL;i++){
        float is_zero=0.;
        for (int j=0;j<FRAME_BLOCK_LEN;j++){
            in_unsort[i][j]=in[INPUT_CHANNEL*j+i];
            is_zero+=in[INPUT_CHANNEL*j+i];
            //printf("%f\t",in[INPUT_CHANNEL*j+i]);
        }
        //printf("\n\n");
        if(is_zero==0.){
            //printf("WHAT THE HECK DID YOU DOl?");
            zero_count++;
            continue;
        }
        //printf("\n---DONE!---\n");
        if(zero_count>=2){
            for(int j=0;j<FRAME_BLOCK_LEN;j++){
                in_sort[ch_tag][j] = in_unsort[i][j];
                //printf("%f\t",in_sort[ch_tag][j]);
                REAL(fft_data[ch_tag],j) = in_sort[ch_tag][j];
                IMAG(fft_data[ch_tag],j) = 0.0;
            }
            //printf("\n\n");
            ch_tag++;
        }
    }
    for(int i=ch_tag;i<USED_CH;i++){
        for(int j=0;j<FRAME_BLOCK_LEN;j++){
            in_sort[i][j] = in_unsort[i-ch_tag][j];
            //printf("%f\t",in_sort[ch_tag][j]);
            REAL(fft_data[i],j) = in_sort[i][j];
            IMAG(fft_data[i],j) = 0.0;
            gsl_matrix_complex_set (in_sort_mat, i, j, double2complex(in_sort[i][j], 0.));
        }
        //printf("\n\n");
    }

    //for(int i=0;i<USED_CH;i++){
    //    for(int j=0;j<FRAME_BLOCK_LEN;j++){
    //        printf("%f\t",in_sort[i][j]);
    //    }
    //    printf("\n\n");
    //}

    for(int i=0;i<USED_CH;i++){
        //printf("\n");

        gsl_fft_complex_radix2_forward(fft_data[i], 1, FRAME_BLOCK_LEN);
        for(int k=0;k<FRAME_BLOCK_LEN;k++){
            gsl_matrix_complex_set (fft_data_mat, i, k, double2complex(REAL(fft_data[i],k), IMAG(fft_data[i],k)));
        }
    }

    //-------Processing in f domain-------
    gsl_matrix_complex *Rxx=gsl_matrix_complex_alloc(fft_data_mat->size1, fft_data_mat->size1);
    gsl_matrix_complex *invR=gsl_matrix_complex_alloc(fft_data_mat->size1, fft_data_mat->size1);
    gsl_matrix_complex *AA=gsl_matrix_complex_alloc(1,USED_CH);
    gsl_matrix_complex *tempww=gsl_matrix_complex_alloc(1,USED_CH);
    gsl_matrix_complex *ww=gsl_matrix_complex_alloc(1,1);
    gsl_matrix_complex *w=gsl_matrix_complex_alloc(USED_CH,1);
    gsl_matrix_complex *y=gsl_matrix_complex_alloc(1, fft_data_mat->size2);
    float doa=0.;
    autocorr(Rxx,fft_data_mat);
    invMat(Rxx,invR);
    for(int i=0;i<USED_CH;i++){
        float tempaa = sin(M_PI/2)*cos(DEG2RAD(doa)-2*i*M_PI)/USED_CH;
        gsl_matrix_complex_set (AA, 0, i, eulers_formula(-1*2*M_PI*Radii*tempaa/LDA));
    }
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, double2complex(1.,0.), AA, invR, double2complex(0.,0.), tempww);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, double2complex(1.,0.), tempww, AA, double2complex(0.,0.), ww);
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_matrix_complex_get(ww,0,0), invR, AA, double2complex(0.,0.), w);
    
    for(int i=0;i<USED_CH;i++){
        gsl_matrix_complex_set(w,i,0,
            gsl_complex_add(
                gsl_complex_mul(
                    gsl_matrix_complex_get(w,i,0),
                    double2complex(Ksnr,0.)),
                double2complex(GAIN*(1-Ksnr*snr), 0.)));
    }
    
    gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, double2complex(1.,0.), w, in_sort_mat, double2complex(0.,0.), y);
    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        REAL(y_arr_f_to_t_domain,i)=gsl_matrix_complex_get(y,0,i).dat[0];
        IMAG(y_arr_f_to_t_domain,i)=gsl_matrix_complex_get(y,0,i).dat[1];
    }
    double original[2*FRAME_BLOCK_LEN];  // original audio using complex form
    double *p_original=original;

    for(int i=0;i<2*FRAME_BLOCK_LEN;i+=2){
        float tmp=0;
        for(int j=0;j<USED_CH;j++){
            tmp+=in_sort[j][i];
        }

        original[i]=tmp;
        original[i+1]=0.;
        *out++ = *y_arr_p++;
        *y_arr_p++;
    }
    gsl_fft_complex_radix2_forward(y_arr_f_to_t_domain, 1, FRAME_BLOCK_LEN);
    snr = SNR(y_arr_f_to_t_domain, 343./LDA, 60.);
    
    //gsl_fft_complex_radix2_inverse(y_arr_f_to_t_domain, 1, FRAME_BLOCK_LEN);

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
