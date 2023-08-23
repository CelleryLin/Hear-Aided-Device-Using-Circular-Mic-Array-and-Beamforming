/*
This program monitors audio in real-time.
Compile and link:

    gcc callback_in_C_wide_Re.c -lportaudio -lgsl -lgslcblas -lfftw3f -lm
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

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define DEG2RAD(x) (x*M_PI/180)

#define FRAME_BLOCK_LEN 512
#define SAMPLING_RATE 44100
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1
#define BANDCOUNT 10

#define TEMPERATURE 340.
#define LDA TEMPERATURE/200.
#define Radii 0.0463

float snr=1.;
float snr0=0.;
float snr_Intergral_val=0.;
float snr_Previous_error=0.;

float doa=1.;
float doa0=0.;
float doa_Intergral_val=0.;
float doa_Previous_error=0.;


PaStream *audioStream;
    
gsl_complex double2complex(float re, float im){
    gsl_complex b;
    b.dat[0]=re;
    b.dat[1]=im;
    return b;
}

float EMA(float a, float a0, float alpha){
    float new_val = a0 + alpha * (a - a0);
    return new_val;
}


float SNR(fftwf_complex a[],double m_freq, double bandwidth){
    // 600, 5.
    float noise_avg=0.;
    float signal_avg=0.;
    float power_spectrum[FRAME_BLOCK_LEN/2+1];
    int m_n = (int)(m_freq*FRAME_BLOCK_LEN/SAMPLING_RATE/2+1);
    int bandwidth_n = (int)(bandwidth*FRAME_BLOCK_LEN/SAMPLING_RATE);

    for (int i=0; i<FRAME_BLOCK_LEN/2+1; i++) {
        power_spectrum[i] = a[i][0] * a[i][0] + a[i][1] * a[i][1];
    }
    
    
    for(int i=0;i<FRAME_BLOCK_LEN/2+1;i++){
        noise_avg+=power_spectrum[i];
    }
    noise_avg/=(FRAME_BLOCK_LEN/2+1);

    for(int i=m_n-bandwidth_n;i<=m_n+bandwidth_n;i++){
        signal_avg+=power_spectrum[i];
    }
    //printf("%f\n", signal_avg);
    signal_avg/=(bandwidth_n*2+1);

    //if((signal_avg/noise_avg)>3){
    //    printf("%f\n",(signal_avg/noise_avg));
    //}
    
    //printf("%f\t%f\n",signal_avg,noise_avg);
    //printf("%f\n",signal_avg/noise_avg);
    return (signal_avg/noise_avg);
    
    //return (signal_avg);
}

void PrintMat(gsl_matrix_complex *m, char info[]){
    printf("Print Matrix %s: \n",info);
    for(int i=0; i<m->size1;i++){
        for(int j=0;j<m->size2;j++){
            printf("%.6f+%.2fi\t", gsl_matrix_complex_get(m,i,j).dat[0],gsl_matrix_complex_get(m,i,j).dat[1]);
        }
        printf("\n");
    }
}

void autocorr(gsl_matrix *Rxx, gsl_matrix *m){
    gsl_matrix *mH=gsl_matrix_alloc(m->size2, m->size1);
    gsl_matrix *Imat=gsl_matrix_alloc(m->size1, m->size1);
    gsl_matrix_set_identity(Imat);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., m, Imat, 0., mH);
    //gsl_matrix_transpose_memcpy(mH, m);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1./m->size2, m, mH,
                    0., Rxx);

    gsl_matrix_free(mH);
    gsl_matrix_free(Imat);

}

void invMat(gsl_matrix *matrix,gsl_matrix *inv){
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);
    int s;
    gsl_linalg_LU_decomp(matrix, p, &s);
    gsl_linalg_LU_invert(matrix, p, inv);
    gsl_permutation_free(p);
}

gsl_complex eulers_formula(double x){    //e^(xi) = cos(x) + i*sin(x)
    gsl_complex a;
    GSL_SET_COMPLEX(&a, cos(x), sin(x));
    return a;
}
float Compress(float n){
    float a, b, tmp;
    a=Ka*n;
    b=Kb*n+1;
    return a/b;
}
void Compress_Mat(gsl_matrix *matrix){
    double a, b, tmp;
    
    for(int i=0; i<matrix->size1; i++){
        for(int j=0; j<matrix->size2; j++){
            if (gsl_matrix_get(matrix, i, j)<=Ka){
                // do nothing
            }
            else{
                gsl_matrix_set(matrix, i, j, 1.);
                //a = gsl_complex_mul(double2complex(Ka, 0.), gsl_matrix_complex_get(matrix, i, j));
                //b = gsl_complex_add(
                //gsl_complex_mul(
                //    double2complex(Kb, 0.),
                //    gsl_matrix_complex_get(matrix, i, j)),
                //double2complex(1., 0.));
                //tmp=gsl_complex_div(a, b);
                //gsl_matrix_complex_set(matrix, i, j, tmp);
            }
        }
    }
}

float hanning(float n, float M){
    return 0.5-0.5*(cos(2*M_PI*n/(M-1)));
}


int audio_callback(const void *inputBuffer,
                   void *outputBuffer, 
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo* timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData){
    
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    float y_arr_f_to_t_domain[FRAME_BLOCK_LEN];
    float *snr_p = &snr;
    float *snr0_p = &snr0;
    float *snr_Intergral_val_p = &snr_Intergral_val;
    float *snr_Previous_error_p = &snr_Previous_error;

    float *doa_p=&doa;
    float *doa0_p=&doa0;
    float *doa_Intergral_val_p=&doa_Intergral_val;
    float *doa_Previous_error_p=&doa_Previous_error;

    static double phase = 0;
    float in_unsort[INPUT_CHANNEL][FRAME_BLOCK_LEN];
    float in_sort[USED_CH][FRAME_BLOCK_LEN];
    fftwf_complex input_fft_data[USED_CH][FRAME_BLOCK_LEN/2+1];
    fftwf_plan plan;
    fftwf_complex output_fft_data[FRAME_BLOCK_LEN/2+1];
    float output[USED_CH][FRAME_BLOCK_LEN];
    int printmat=0;
    gsl_matrix_complex *fft_data_mat=gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN/2+1);
    gsl_matrix *in_sort_mat=gsl_matrix_alloc(USED_CH, FRAME_BLOCK_LEN);
    gsl_matrix_complex *insort_mat_complex = gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    //printf("\n\n-----------start print------------\n");
    int ch_tag=0;
    int zero_count=0;
    for (int i=0;i<INPUT_CHANNEL;i++){
        float is_zero=0.;
        for (int j=0;j<FRAME_BLOCK_LEN;j++){
            in_unsort[i][j]=in[INPUT_CHANNEL*j+i]; //*hanning(j, FRAME_BLOCK_LEN);
            is_zero+=in[INPUT_CHANNEL*j+i];
            //printf("%f\t",in[INPUT_CHANNEL*j+i]);
        }
        //printf("\n\n");
        if(is_zero==0. && zero_count!=2){
            //printf("WHAT THE HECK DID YOU DOl?");
            zero_count++;
            continue;
        }
        //printf("\n---DONE!---\n");
        if(zero_count>=2){
            for(int j=0;j<FRAME_BLOCK_LEN;j++){
                in_sort[ch_tag][j] = in_unsort[i][j];
            }
            //printf("\n\n");
            ch_tag++;
        }
    }
    for(int i=ch_tag;i<USED_CH;i++){
        for(int j=0;j<FRAME_BLOCK_LEN;j++){
            in_sort[i][j] = in_unsort[i-ch_tag][j];
        }
        //printf("\n\n");
    }
    for(int i=0;i<USED_CH;i++){
        for(int j=0;j<FRAME_BLOCK_LEN;j++){
            gsl_matrix_set (in_sort_mat, i, j, (double)in_sort[i][j]);
            gsl_matrix_complex_set (insort_mat_complex, i, j, double2complex(in_sort[i][j], 0.));
        }
        //printf("\n\n");
    }

    // for(int i=0;i<USED_CH;i++){
    //     for(int j=0;j<FRAME_BLOCK_LEN;j++){
    //         printf("%f\t",in_sort[i][j]);
    //     }
    //     printf("\n\n");
    // }

    for(int i=0;i<USED_CH;i++){
        //printf("\n");
        plan = fftwf_plan_dft_r2c_1d(FRAME_BLOCK_LEN, in_sort[i], input_fft_data[i], FFTW_ESTIMATE);
        fftwf_execute(plan);


        //gsl_fft_complex_radix2_forward(input_fft_data[i], 1, FRAME_BLOCK_LEN);
        for(int k=0;k<FRAME_BLOCK_LEN/2+1;k++){
            gsl_matrix_complex_set (fft_data_mat, i, k, double2complex(input_fft_data[i][k][0], input_fft_data[i][k][1]));
        }
    }
    
    //for(int i=0;i<FRAME_BLOCK_LEN;i++){
    //    printf("%f\t",input_fft_data[4][i][0]);
    //}
    //printf("\n");
    
    //PrintMat(fft_data_mat,"---");
    //-------Processing in f domain-------
    gsl_matrix_complex *invR=gsl_matrix_complex_alloc(fft_data_mat->size1, fft_data_mat->size1);


    gsl_matrix *Rxx_Re=gsl_matrix_alloc(in_sort_mat->size1, in_sort_mat->size1);
    gsl_matrix *invR_Re=gsl_matrix_alloc(in_sort_mat->size1, in_sort_mat->size1);


    gsl_matrix_complex *AA=gsl_matrix_complex_alloc(1,USED_CH);
    gsl_matrix_complex *bestAA=gsl_matrix_complex_alloc(1,USED_CH);
    gsl_matrix_complex *tempww=gsl_matrix_complex_alloc(1,USED_CH);
    gsl_matrix_complex *ww=gsl_matrix_complex_alloc(1,1);
    gsl_matrix_complex *bestww=gsl_matrix_complex_alloc(1,1);
    gsl_matrix_complex *w=gsl_matrix_complex_alloc(USED_CH,1);
    gsl_matrix_complex *y=gsl_matrix_complex_alloc(1, FRAME_BLOCK_LEN);
    gsl_matrix_complex_set_zero(y);
    double max_ww=-1;
    autocorr(Rxx_Re,in_sort_mat);
    invMat(Rxx_Re,invR_Re);
    //Compress_Mat(invR_Re);
    //PrintMat(invR,"INVR");
    //abort();
    
    for(int idoa=0; idoa<=360; idoa+=45){
        for(int i=0;i<USED_CH;i++){
            float tempaa = (float) (sin(M_PI/2)*cos(DEG2RAD((double)idoa)-2*i*M_PI/USED_CH));
            gsl_matrix_complex_set (AA, 0, i, eulers_formula(-1*2*M_PI*Radii*tempaa/(LDA)));
            for(int j=0;j<USED_CH;j++){
                gsl_matrix_complex_set (invR, i, j, double2complex(gsl_matrix_get(invR_Re, i, j), 0.));
            }
        }
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, double2complex(1.,0.), AA, invR, double2complex(0.,0.), tempww);
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, double2complex(1.,0.), tempww, AA, double2complex(0.,0.), ww);
        //gsl_matrix_complex_set(ww,0,0,double2complex(
        //        Compress(gsl_matrix_complex_get(ww,0,0).dat[0]),    //Compress the gain of ww
        //        gsl_matrix_complex_get(ww,0,0).dat[1]));
        //printf("%f\t%f\n",gsl_complex_div(double2complex(1.,0.),gsl_matrix_complex_get(ww,0,0)));
        //PrintMat(ww,"ww");
        if(gsl_matrix_complex_get(ww,0,0).dat[0]>max_ww){
            *doa_p=idoa;
            max_ww=gsl_matrix_complex_get(ww,0,0).dat[0];
            gsl_matrix_complex_memcpy(bestww,ww);
            //gsl_matrix_complex_memcpy(bestAA,AA);
        }
    }

    //*doa_p=PID(*doa_p,*doa0_p,0.1,0,0,*doa_Intergral_val_p,*doa_Previous_error_p);
    *doa_p = EMA(*doa_p, *doa0_p, 0.1);
    *doa0_p=*doa_p;
    // *doa_p=0;
    // for(int i=0;i<*doa_p/5;i++){
    //     printf(" ");
    // }
    // printf("|\n");


    for(int l=start_freq; l<end_freq; l+=100){    //Wide band
        for(int i=0;i<USED_CH;i++){
                float tempaa = (float) (sin(M_PI/2)*cos(DEG2RAD((double)*doa_p)-2*i*M_PI/USED_CH));
                gsl_matrix_complex_set (bestAA, 0, i, eulers_formula(-1*2*M_PI*Radii*tempaa/(TEMPERATURE/l)));
        }

        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_div(double2complex(1.,0.),gsl_matrix_complex_get(bestww,0,0)), invR, bestAA, double2complex(0.,0.), w);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, double2complex(GAIN*(*snr_p)/(float)BANDCOUNT,0.), w, insort_mat_complex, double2complex(1./(float)BANDCOUNT,0.), y);
        //gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, double2complex(GAIN_BIG/(float)BANDCOUNT,0.), w, insort_mat_complex, double2complex(1./(float)BANDCOUNT,0.), y);
    }
    //PrintMat(y,"y");
    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        //y_arr_f_to_t_domain[i] = (float) gsl_complex_abs(gsl_matrix_complex_get(y,0,i));
        y_arr_f_to_t_domain[i] = (float) gsl_matrix_complex_get(y,0,i).dat[0];
        //printf("%f\t",y_arr_f_to_t_domain[i]);
    }
    //printf("%f\t",y_arr_f_to_t_domain[0]);
    //printf("%f\n",*snr_p);
    float original[FRAME_BLOCK_LEN];  // original audio using complex form
    float *p_original=original;
    float *y_arr_p = y_arr_f_to_t_domain;
    float y_arr_f_to_t_domain_filtered[FRAME_BLOCK_LEN];
    float *y_arr_filtered_p = y_arr_f_to_t_domain_filtered;


    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        // Savitzky–Golay filter Y = sum(Cy), C = [-3, 12, 17, 12, -3]/35
        y_arr_f_to_t_domain_filtered[i] = 0;
        if(i-3>=0){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i-3]*-2.;
        }
        if(i-2>=0){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i-2]*3.;
        }
        if(i-1>=0){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i-1]*6.;
        }
        y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i]*7.;
        if(i+1<FRAME_BLOCK_LEN){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i+1]*6.;
        }
        if(i+2<FRAME_BLOCK_LEN){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i+2]*3.;
        }
        if(i+3<FRAME_BLOCK_LEN){
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i+3]*-2.;
        }
        y_arr_f_to_t_domain_filtered[i] /= 21.;
        // float tmp=0;
        // for(int j=0;j<USED_CH;j++){
        //     tmp+=(in_sort[j][i]/6);
        //     //tmp+=output[j][i];
        // }
        // original[i]=tmp;
        //*out++ = *p_original++;
        *out++ = *y_arr_filtered_p++;
        //*out++ = *y_arr_p++;
        //*y_arr_p++;
    }
    //rnnoise_process_frame(st, out, out);
    plan = fftwf_plan_dft_r2c_1d(FRAME_BLOCK_LEN, y_arr_f_to_t_domain, output_fft_data, FFTW_ESTIMATE);
    fftwf_execute(plan);
    *snr_p = SNR(output_fft_data, 600, 10.);
    *snr_p = EMA(*snr_p,*snr0_p, 0.01);
    if (*snr_p<=snr_gate) *snr_p=min_snr;
    //printf("%f\n",*snr_p);
    *snr0_p = *snr_p;
    
    // printf("%f\n",*snr_p);
    //gsl_fft_complex_radix2_inverse(y_arr_f_to_t_domain, 1, FRAME_BLOCK_LEN);

    //Free these all Fucking things
    gsl_matrix_complex_free(invR);
    gsl_matrix_free(Rxx_Re);
    gsl_matrix_free(invR_Re);
    gsl_matrix_complex_free(AA);
    gsl_matrix_complex_free(tempww);
    gsl_matrix_complex_free(ww);
    gsl_matrix_complex_free(w);
    gsl_matrix_complex_free(y);
    gsl_matrix_complex_free(fft_data_mat);
    gsl_matrix_free(in_sort_mat);
    gsl_matrix_complex_free(insort_mat_complex);
    gsl_matrix_complex_free(bestww);
    gsl_matrix_complex_free(bestAA);
    fftwf_destroy_plan(plan);

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
