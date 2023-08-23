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

int audio_callback(const void *inputBuffer,
                   void *outputBuffer,
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData)
{

    float *in = (float *)inputBuffer, *out = (float *)outputBuffer;
    float y_arr_f_to_t_domain[FRAME_BLOCK_LEN];
    float *snr_p = &snr;
    float *snr0_p = &snr0;
    float *snr_Intergral_val_p = &snr_Intergral_val;
    float *snr_Previous_error_p = &snr_Previous_error;

    float *doa_p = &doa;
    float *doa0_p = &doa0;

    static double phase = 0;
    float in_unsort[INPUT_CHANNEL][FRAME_BLOCK_LEN];
    float in_sort[USED_CH][FRAME_BLOCK_LEN];
    float in_sort_window[USED_CH][FRAME_BLOCK_LEN];
    fftwf_complex input_fft_data[USED_CH][FRAME_BLOCK_LEN / 2 + 1];
    fftwf_plan plan;
    fftwf_complex output_fft_data[FRAME_BLOCK_LEN / 2 + 1];
    float output[USED_CH][FRAME_BLOCK_LEN];
    int printmat = 0;
    gsl_matrix_complex *fft_data_mat = gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN / 2 + 1);
    gsl_matrix_complex *in_sort_mat = gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    gsl_matrix_complex *in_sort_window_mat = gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    // printf("\n\n-----------start print------------\n");
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
            // printf("WHAT THE HECK DID YOU DO?");
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
    for (int i = 0; i < USED_CH; i++)
    {
        for (int j = 0; j < FRAME_BLOCK_LEN; j++)
        {
            gsl_matrix_complex_set(in_sort_mat, i, j, double2complex(in_sort[i][j], 0.));
            // gsl_matrix_complex_set(in_sort_window_mat, i, j, double2complex(in_sort_window[i][j], 0.));
        }
        // printf("\n\n");
    }

    // for(int i=0;i<USED_CH;i++){
    //     for(int j=0;j<FRAME_BLOCK_LEN;j++){
    //         printf("%f\t",in_sort[i][j]);
    //     }
    //     printf("\n\n");
    // }

    for (int i = 0; i < USED_CH; i++)
    {
        // printf("\n");
        plan = fftwf_plan_dft_r2c_1d(FRAME_BLOCK_LEN, in_sort[i], input_fft_data[i], FFTW_ESTIMATE);
        fftwf_execute(plan);

        // gsl_fft_complex_radix2_forward(input_fft_data[i], 1, FRAME_BLOCK_LEN);
        for (int k = 0; k < FRAME_BLOCK_LEN / 2 + 1; k++)
        {
            gsl_matrix_complex_set(fft_data_mat, i, k, double2complex(input_fft_data[i][k][0], input_fft_data[i][k][1]));
        }
    }
    float vad = 0.;
    vad = energy_vad(input_fft_data[0], vad_threshold);
    // printf("vad: %f\n", vad);

    // for(int i=0;i<FRAME_BLOCK_LEN;i++){
    //     printf("%f\t",input_fft_data[4][i][0]);
    // }
    // printf("\n");

    // PrintMat(fft_data_mat,"---");
    //-------Processing in f domain-------
    gsl_matrix_complex *Rxx = gsl_matrix_complex_alloc(fft_data_mat->size1, fft_data_mat->size1);
    gsl_matrix_complex *invR = gsl_matrix_complex_alloc(fft_data_mat->size1, fft_data_mat->size1);
    gsl_matrix_complex *AA = gsl_matrix_complex_alloc(1, USED_CH);
    gsl_matrix_complex *bestAA = gsl_matrix_complex_alloc(1, USED_CH);
    gsl_matrix_complex *tempww = gsl_matrix_complex_alloc(1, USED_CH);
    gsl_matrix_complex *ww = gsl_matrix_complex_alloc(1, 1);
    gsl_matrix_complex *bestww = gsl_matrix_complex_alloc(1, 1);
    gsl_matrix_complex *w = gsl_matrix_complex_alloc(USED_CH, 1);
    gsl_matrix_complex *y = gsl_matrix_complex_alloc(1, FRAME_BLOCK_LEN);
    gsl_matrix_complex_set_zero(y);
    double max_ww = -1;
    autocorr(Rxx, fft_data_mat);
    invMat(Rxx, invR);
    // Compress_Mat(invR);
    // PrintMat(invR,"INVR");
    // abort();

    for (int idoa = 0; idoa <= 360; idoa += 45)
    {
        for (int i = 0; i < USED_CH; i++)
        {
            float tempaa = (float)(sin(M_PI / 2) * cos(DEG2RAD((double)idoa) - 2 * i * M_PI / USED_CH));
            gsl_matrix_complex_set(AA, 0, i, eulers_formula(-1 * 2 * M_PI * Radii * tempaa / (LDA)));
        }
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, double2complex(1., 0.), AA, invR, double2complex(0., 0.), tempww);
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, double2complex(1., 0.), tempww, AA, double2complex(0., 0.), ww);
        // gsl_matrix_complex_set(ww,0,0,double2complex(
        //         Compress(gsl_matrix_complex_get(ww,0,0).dat[0]),    //Compress the gain of ww
        //         gsl_matrix_complex_get(ww,0,0).dat[1]));
        // printf("%f\t%f\n",gsl_complex_div(double2complex(1.,0.),gsl_matrix_complex_get(ww,0,0)));
        // PrintMat(ww,"ww");
        if (gsl_matrix_complex_get(ww, 0, 0).dat[0] > max_ww)
        {
            *doa_p = idoa;
            max_ww = gsl_matrix_complex_get(ww, 0, 0).dat[0];
            gsl_matrix_complex_memcpy(bestww, ww);
            // gsl_matrix_complex_memcpy(bestAA,AA);
        }
    }

    //*doa_p=PID(*doa_p,*doa0_p,0.1,0,0,*doa_Intergral_val_p,*doa_Previous_error_p);
    *doa_p = EMA(*doa_p, *doa0_p, EMA_RATE);
    *doa0_p = *doa_p;
    // for(int i=0;i<*doa_p/5;i++){
    //    printf(" ");
    //}
    // printf("|\n");

    for (int l = start_freq; l < end_freq; l += 100)
    { // Wide band
        for (int i = 0; i < USED_CH; i++)
        {
            float tempaa = (float)(sin(M_PI / 2) * cos(DEG2RAD((double)*doa_p) - 2 * i * M_PI / USED_CH));
            gsl_matrix_complex_set(bestAA, 0, i, eulers_formula(-1 * 2 * M_PI * Radii * tempaa / (TEMPERATURE / l)));
        }

        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_div(double2complex(1., 0.), gsl_matrix_complex_get(bestww, 0, 0)), invR, bestAA, double2complex(0., 0.), w);
        //gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, double2complex(GAIN * (*snr_p) * vad / (float)BANDCOUNT, 0.), w, in_sort_mat, double2complex(1. / (float)BANDCOUNT, 0.), y);
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, double2complex(GAIN_BIG*vad/(float)BANDCOUNT,0.), w, in_sort_mat, double2complex(1./(float)BANDCOUNT,0.), y);
    }
    // PrintMat(y,"y");
    for (int i = 0; i < FRAME_BLOCK_LEN; i++)
    {
        // y_arr_f_to_t_domain[i] = (float) gsl_complex_abs(gsl_matrix_complex_get(y,0,i));
        y_arr_f_to_t_domain[i] = (float)gsl_matrix_complex_get(y, 0, i).dat[0];
        // printf("%f\t",y_arr_f_to_t_domain[i]);
    }
    // printf("%f\t",y_arr_f_to_t_domain[0]);
    // printf("%f\n",*snr_p);
    float original[FRAME_BLOCK_LEN]; // original audio using complex form
    float *p_original = original;
    float *y_arr_p = y_arr_f_to_t_domain;
    float y_arr_f_to_t_domain_filtered[FRAME_BLOCK_LEN];
    float *y_arr_filtered_p = y_arr_f_to_t_domain_filtered;

    for (int i = 0; i < FRAME_BLOCK_LEN; i++)
    {
        // Savitzky–Golay filter Y = sum(Cy), C = [-3, 12, 17, 12, -3]/35
        y_arr_f_to_t_domain_filtered[i] = 0;
        if (i - 3 >= 0)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i - 3] * -2.;
        }
        if (i - 2 >= 0)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i - 2] * 3.;
        }
        if (i - 1 >= 0)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i - 1] * 6.;
        }
        y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i] * 7.;
        if (i + 1 < FRAME_BLOCK_LEN)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i + 1] * 6.;
        }
        if (i + 2 < FRAME_BLOCK_LEN)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i + 2] * 3.;
        }
        if (i + 3 < FRAME_BLOCK_LEN)
        {
            y_arr_f_to_t_domain_filtered[i] += y_arr_f_to_t_domain[i + 3] * -2.;
        }
        y_arr_f_to_t_domain_filtered[i] /= 21.;
        // float tmp=0;
        // for(int j=0;j<USED_CH;j++){
        //     tmp+=(in_sort[j][i]/6);
        //     //tmp+=output[j][i];
        // }
        // original[i]=tmp;
        // *out++ = *p_original++;
        *out++ = *y_arr_filtered_p++;
        // *out++ = *y_arr_p++;
        //*y_arr_p++;
    }
    // rnnoise_process_frame(st, out, out);
    plan = fftwf_plan_dft_r2c_1d(FRAME_BLOCK_LEN, y_arr_f_to_t_domain, output_fft_data, FFTW_ESTIMATE);
    fftwf_execute(plan);
    //*snr_p = SNR(input_fft_data[0], snr_start_freq, snr_end_freq);
    *snr_p = SNR(output_fft_data, snr_start_freq, snr_end_freq);
    // printf("%f\n",*snr_p);
    //  *snr_p = PID(*snr_p,*snr0_p, 0.001, 0, 0, *snr_Intergral_val_p, *snr_Previous_error_p);
    *snr_p = EMA(*snr_p, *snr0_p, 0.1);
    if (*snr_p < min_snr)
    {
        *snr_p = min_snr;
    }

    *snr0_p = *snr_p;

    // printf("%f\n",*snr_p);
    // gsl_fft_complex_radix2_inverse(y_arr_f_to_t_domain, 1, FRAME_BLOCK_LEN);

    // Free these all Fucking things
    gsl_matrix_complex_free(Rxx);
    gsl_matrix_complex_free(invR);
    gsl_matrix_complex_free(AA);
    gsl_matrix_complex_free(tempww);
    gsl_matrix_complex_free(ww);
    gsl_matrix_complex_free(w);
    gsl_matrix_complex_free(y);
    gsl_matrix_complex_free(fft_data_mat);
    gsl_matrix_complex_free(in_sort_mat);
    gsl_matrix_complex_free(bestww);
    gsl_matrix_complex_free(bestAA);
    fftwf_destroy_plan(plan);

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

void terminate_stuff()
{

    Pa_StopStream(audioStream);  /* stop the callback mechanism */
    Pa_CloseStream(audioStream); /* destroy the audio stream object */
    Pa_Terminate();              /* terminate portaudio */
    only_bf_terminate();         /* terminate beamforming unit */
}

int main()
{
    init_stuff();
    while (getchar() != ' ')
        Pa_Sleep(100);
    terminate_stuff();
    return 0;
}
