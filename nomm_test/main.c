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

#include "portaudio.h"

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define DEG2RAD(x) (x*M_PI/180)

#define FRAME_BLOCK_LEN 1024
#define SAMPLING_RATE 44100
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1
#define BANDCOUNT 10

#define TEMPERATURE 340.
#define LDA TEMPERATURE/200.
#define Radii 0.0463

#define Ka 100000.           //Adjust compressor strength
#define Kb 1.          //Adjust compressor strength (Bypass: Ka=1, Kb=0)
#define Kmvdr 1.        //Adjust snr strength
#define snr_gate 0          //treshold of snr
#define GAIN 1

#define wKp 0.7
#define wKi 0.0
#define wKd 0.

// --- rnn denoise ---
#include "nnom.h"
#include "denoise_weights.h"
#include "mfcc.h"

 // the bandpass filter coefficiences
#include "equalizer_coeff.h" 

#define NUM_FEATURES NUM_FILTER

#define _MAX(x, y) (((x) > (y)) ? (x) : (y))
#define _MIN(x, y) (((x) < (y)) ? (x) : (y))
#define SaturaLH(N, L, H) (((N)<(L))?(L):(((N)>(H))?(H):(N)))

#define RNN_GAIN_IN     5.0f
#define RNN_GAIN_OUT    0.5f
#define NUM_CHANNELS 	1
#define SAMPLE_RATE 	SAMPLING_RATE
#define AUDIO_FRAME_LEN FRAME_BLOCK_LEN

// audio buffer for input
float audio_buffer[AUDIO_FRAME_LEN] = {0};
int16_t audio_buffer_16bit[AUDIO_FRAME_LEN] = {0};

// buffer for output
int16_t audio_buffer_filtered[AUDIO_FRAME_LEN] = { 0 };

// mfcc features and their derivatives
float mfcc_feature[NUM_FEATURES] = { 0 };
float mfcc_feature_prev[NUM_FEATURES] = { 0 };
float mfcc_feature_diff[NUM_FEATURES] = { 0 };
float mfcc_feature_diff_prev[NUM_FEATURES] = { 0 };
float mfcc_feature_diff1[NUM_FEATURES] = { 0 };

// features for NN
float nn_features[64] = {0};
int8_t nn_features_q7[64] = {0};

// NN results, which is the gains for each frequency band
float band_gains[NUM_FILTER] = {0};
float band_gains_prev[NUM_FILTER] = {0};

// 0db gains coefficient
float coeff_b[NUM_FILTER][NUM_COEFF_PAIR] = FILTER_COEFF_B;
float coeff_a[NUM_FILTER][NUM_COEFF_PAIR] = FILTER_COEFF_A;

// dynamic gains coefficient
float b_[NUM_FILTER][NUM_COEFF_PAIR] = {0};

// nnom model
nnom_model_t *model;
mfcc_t *mfcc;

// -------------------

float snr=1.;
float snr0=0.;
float snr_Intergral_val=0.;
float snr_Previous_error=0.;

float doa=1.;
float doa0=0.;
float doa_Intergral_val=0.;
float doa_Previous_error=0.;

gsl_matrix_complex *w0;
gsl_matrix_complex *w_Intergral_val;
gsl_matrix_complex *w_Previous_error;

PaStream *audioStream;
    
gsl_complex double2complex(float re, float im){
    gsl_complex b;
    b.dat[0]=re;
    b.dat[1]=im;
    return b;
}

gsl_matrix_complex PID_matrix(gsl_matrix_complex *a, gsl_matrix_complex *a0){
    gsl_matrix_complex *Derivative_val=gsl_matrix_complex_alloc(a->size1, a->size2);
    gsl_matrix_complex *error=gsl_matrix_complex_alloc(a->size1, a->size2);
    gsl_matrix_complex *new_val=gsl_matrix_complex_alloc(a->size1, a->size2);
    for(int i=0;i<a->size1;i++){
        for(int j=0;j<a->size2;j++){
            gsl_matrix_complex_set(error,i,j,gsl_complex_sub(
                gsl_matrix_complex_get(a,i,j),
                gsl_matrix_complex_get(a0,i,j)));
            gsl_matrix_complex_set(w_Intergral_val,i,j,gsl_complex_add(
                gsl_matrix_complex_get(w_Intergral_val,i,j),
                gsl_matrix_complex_get(error,i,j)));
            gsl_matrix_complex_set(Derivative_val,i,j,gsl_complex_sub(
                gsl_matrix_complex_get(error,i,j),
                gsl_matrix_complex_get(w_Previous_error,i,j)));
            gsl_matrix_complex_set(w_Previous_error,i,j,gsl_matrix_complex_get(error,i,j));
            
            gsl_matrix_complex_set(new_val,i,j,
                gsl_complex_add(gsl_complex_add(gsl_complex_add(
                    gsl_matrix_complex_get(a0,i,j),
                    gsl_complex_mul(double2complex(wKp,0.), gsl_matrix_complex_get(error,i,j))),
                    gsl_complex_mul(double2complex(wKi,0.), gsl_matrix_complex_get(w_Intergral_val,i,j))),
                    gsl_complex_mul(double2complex(wKd,0.), gsl_matrix_complex_get(Derivative_val,i,j))));
        }
    }
    gsl_matrix_complex_free(Derivative_val);
    gsl_matrix_complex_free(error);
    return *new_val;
}

float PID(float a, float a0, float Kp, float Ki, float Kd, float Intergral_val, float Previous_error){
    float error=a-a0;
    snr_Intergral_val += error;
    float Derivative_val = error - snr_Previous_error;
    snr_Previous_error = error;
    double new_val =a0 + (Kp * error + Ki * snr_Intergral_val + Kd * Derivative_val);
    return new_val;
}

float EMA(float a, float a0, float alpha){
    float new_val = a0 + alpha * (a - a0);
    return new_val;
}

float SNR(fftwf_complex a[],double m_freq, double bandwidth){
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
    if(signal_avg/noise_avg>=snr_gate){
        return (signal_avg/noise_avg);
    }
    else{
        return (0);
    }
    
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

void autocorr(gsl_matrix_complex *Rxx, gsl_matrix_complex *m){
    gsl_matrix_complex *mH=gsl_matrix_complex_alloc(m->size2, m->size1);
    gsl_matrix_complex *Imat=gsl_matrix_complex_alloc(m->size1, m->size1);
    gsl_matrix_complex_set_identity(Imat);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,double2complex(1.,0.),m,Imat,double2complex(0.,0.),mH);
    //gsl_matrix_transpose_memcpy(mH, m);
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
                    double2complex(1./m->size2,0.), m, mH,
                    double2complex(0.,0.), Rxx);

    gsl_matrix_complex_free(mH);
    gsl_matrix_complex_free(Imat);

}

void invMat(gsl_matrix_complex *matrix,gsl_matrix_complex *inv){
    gsl_permutation *p = gsl_permutation_alloc(matrix->size1);
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
float Compress(float n){
    float a, b, tmp;
    a=Ka*n;
    b=Kb*n+1;
    return a/b;
}
void Compress_Mat(gsl_matrix_complex *matrix){
    gsl_complex a, b, tmp;
    
    for(int i=0; i<matrix->size1; i++){
        for(int j=0; j<matrix->size2; j++){
            if (gsl_matrix_complex_get(matrix, i, j).dat[0]<=Ka){
                // do nothing
            }
            else{
                gsl_matrix_complex_set(matrix, i, j, double2complex(1., 0.));
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

// ---------------------------------------------------------
// -----------------  rnn_denoise_function -----------------
// ---------------------------------------------------------

void y_h_update(float *y_h, uint32_t len)
{
	for (uint32_t i = len-1; i >0 ;i--)
		y_h[i] = y_h[i-1];
}

//  equalizer by multiple n order iir band pass filter. 
// y[i] = b[0] * x[i] + b[1] * x[i - 1] + b[2] * x[i - 2] - a[1] * y[i - 1] - a[2] * y[i - 2]...
void equalizer(float* x, float* y, uint32_t signal_len, float *b, float *a, uint32_t num_band, uint32_t num_order)
{
	// the y history for each band
	static float y_h[NUM_FILTER][NUM_COEFF_PAIR] = { 0 };
	static float x_h[NUM_COEFF_PAIR * 2] = { 0 };
	uint32_t num_coeff = num_order * 2 + 1;

	// i <= num_coeff (where historical x is involved in the first few points)
	// combine state and new data to get a continual x input. 
	memcpy(x_h + num_coeff, x, num_coeff * sizeof(float));
	for (uint32_t i = 0; i < num_coeff; i++)
	{
		y[i] = 0;
		for (uint32_t n = 0; n < num_band; n++)
		{
			y_h_update(y_h[n], num_coeff);
			y_h[n][0] = b[n * num_coeff] * x_h[i+ num_coeff];
			for (uint32_t c = 1; c < num_coeff; c++)
				y_h[n][0] += b[n * num_coeff + c] * x_h[num_coeff + i - c] - a[n * num_coeff + c] * y_h[n][c];
			y[i] += y_h[n][0];
		}
	}
	// store the x for the state of next round
	memcpy(x_h, &x[signal_len - num_coeff], num_coeff * sizeof(float));
	
	// i > num_coeff; the rest data not involed the x history
	for (uint32_t i = num_coeff; i < signal_len; i++)
	{
		y[i] = 0;
		for (uint32_t n = 0; n < num_band; n++)
		{
			y_h_update(y_h[n], num_coeff);
			y_h[n][0] = b[n * num_coeff] * x[i];
			for (uint32_t c = 1; c < num_coeff; c++)
				y_h[n][0] += b[n * num_coeff + c] * x[i - c] - a[n * num_coeff + c] * y_h[n][c];
			y[i] += y_h[n][0];
		}	
	}
}

// set dynamic gains. Multiple gains x b_coeff
void set_gains(float *b_in, float *b_out,  float* gains, uint32_t num_band, uint32_t num_order)
{
	uint32_t num_coeff = num_order * 2 + 1;
	for (uint32_t i = 0; i < num_band; i++)
		for (uint32_t c = 0; c < num_coeff; c++)
			b_out[num_coeff *i + c] = b_in[num_coeff * i + c] * gains[i]; 
}

void quantize_data(float*din, int8_t *dout, uint32_t size, uint32_t int_bit)
{
	float limit = (1 << int_bit); 
	for(uint32_t i=0; i<size; i++)
		dout[i] = (int8_t)(_MAX(_MIN(din[i], limit), -limit) / limit * 127);
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// ----------------------------------------------------------------



int audio_callback(const void *inputBuffer,
                   void *outputBuffer, 
                   unsigned long framesPerBuffer,
                   const PaStreamCallbackTimeInfo* timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData){
    
    float *in = (float*) inputBuffer, *out = (float*)outputBuffer;
    float y_arr_f_to_t_domain[FRAME_BLOCK_LEN];
    float *y_arr_p = y_arr_f_to_t_domain;
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
    fftwf_complex output_fft_data[FRAME_BLOCK_LEN/2+1];
    float output[USED_CH][FRAME_BLOCK_LEN];
    int printmat=0;
    gsl_matrix_complex *fft_data_mat=gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN/2+1);
    gsl_matrix_complex *in_sort_mat=gsl_matrix_complex_alloc(USED_CH, FRAME_BLOCK_LEN);
    int32_t *p_new_data;
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
    float original[FRAME_BLOCK_LEN];  // original audio using complex form
    float *p_original=original;

    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        float tmp=0;
        for(int j=0;j<USED_CH;j++){
            tmp+=(in_sort[j][i]/USED_CH);
            //tmp+=output[j][i];
        }
        original[i]=tmp;
        //*out++ = *p_original++;
        //*out++ = *y_arr_p++;
        //*y_arr_p++;
    }

    // -------------------------------------------
    // --------------- rnn denoise ---------------
    // -------------------------------------------



    // memcpy(audio_buffer_16bit, &audio_buffer_16bit[AUDIO_FRAME_LEN/2], AUDIO_FRAME_LEN/2*sizeof(int16_t));

    for(int i = 0; i < AUDIO_FRAME_LEN; i++){
        audio_buffer_16bit[i] = (int16_t)SaturaLH((original[i] * 32768.f * RNN_GAIN_IN), -32768.f, 32767.f);
    }


    // get mfcc
    mfcc_compute(mfcc, audio_buffer_16bit, mfcc_feature);
    
    // get the first and second derivative of mfcc
    for(uint32_t i=0; i< NUM_FEATURES; i++){
        mfcc_feature_diff[i] = mfcc_feature[i] - mfcc_feature_prev[i];
        mfcc_feature_diff1[i] = mfcc_feature_diff[i] - mfcc_feature_diff_prev[i];
    }
    memcpy(mfcc_feature_prev, mfcc_feature, NUM_FEATURES * sizeof(float));
    memcpy(mfcc_feature_diff_prev, mfcc_feature_diff, NUM_FEATURES * sizeof(float));
    
    // combine MFCC with derivatives for the NN features
    memcpy(nn_features, mfcc_feature, NUM_FEATURES*sizeof(float));
    memcpy(&nn_features[NUM_FEATURES], mfcc_feature_diff, 10*sizeof(float));
    memcpy(&nn_features[NUM_FEATURES+10], mfcc_feature_diff1, 10*sizeof(float));

    // quantise them using the same scale as training data (in keras), by 2^n. 
    quantize_data(nn_features, nn_features_q7, NUM_FEATURES+20, 3);
    
    // run the mode with the new input
    memcpy(nnom_input_data, nn_features_q7, sizeof(nnom_input_data));
    model_run(model);

    // read the result, convert it back to float (q0.7 to float)
    for(int i=0; i< NUM_FEATURES; i++)
        band_gains[i] = (float)(nnom_output_data[i]) / 127.f;
    
    // one more step, limit the change of gians, to smooth the speech, per RNNoise paper
    for(int i=0; i< NUM_FEATURES; i++)
        band_gains[i] = _MAX(band_gains_prev[i]*0.8f, band_gains[i]); 
    memcpy(band_gains_prev, band_gains, NUM_FEATURES *sizeof(float));
    
    // update filter coefficient to applied dynamic gains to each frequency band 
    set_gains((float*)coeff_b, (float*)b_, band_gains, NUM_FILTER, NUM_ORDER);

    // convert 16bit to float for equalizer
    for (int i = 0; i < AUDIO_FRAME_LEN; i++)
        audio_buffer[i] = audio_buffer_16bit[i] / 32768.f;
            
    // finally, we apply the equalizer to this audio frame to denoise
    equalizer(audio_buffer, &audio_buffer[AUDIO_FRAME_LEN], AUDIO_FRAME_LEN, (float*)b_,(float*)coeff_a, NUM_FILTER, NUM_ORDER);
    
    // convert it back to int16
    for (int i = 0; i < AUDIO_FRAME_LEN; i++)
        audio_buffer_filtered[i] = audio_buffer[i] * 32768.f *0.7f; // 0.7 is the filter band overlapping factor


    
    float audio_buffer_filtered_float[FRAME_BLOCK_LEN];
    float *p_audio_buffer_filtered_float=audio_buffer_filtered_float;

    for(int i=0;i<FRAME_BLOCK_LEN;i++){
        audio_buffer_filtered_float[i] = (float)audio_buffer_filtered[i];
        audio_buffer_filtered_float[i] /= 32768.f;
        audio_buffer_filtered_float[i] *= RNN_GAIN_OUT;
        *out++ = *p_audio_buffer_filtered_float++;
        // *out++ = *p_original++;
        //*out++ = *y_arr_p++;
        //*y_arr_p++;
    }
    // if(nnom_output_data1[0] >= 64)
    //     printf("speech\n");
    // else
    //     printf("noise\n");



    // ------------------------------------------------
    // ------------------------------------------------
    // ------------------------------------------------



    gsl_matrix_complex_free(fft_data_mat);
    gsl_matrix_complex_free(in_sort_mat);

    return paContinue;
}

void init_stuff(){
    int i,id;
    const PaDeviceInfo *info;
    const PaHostApiInfo *hostapi;
    PaStreamParameters outputParameters, inputParameters;
    w0 = gsl_matrix_complex_alloc(USED_CH, 1);
    w_Intergral_val = gsl_matrix_complex_alloc(USED_CH, 1);
    w_Previous_error = gsl_matrix_complex_alloc(USED_CH, 1);
    gsl_matrix_complex_set_zero(w0);
    gsl_matrix_complex_set_zero(w_Intergral_val);
    gsl_matrix_complex_set_zero(w_Previous_error);

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
    printf("load model\n");
    model = nnom_model_create();
    mfcc = mfcc_create(NUM_FEATURES, 0, NUM_FEATURES, 512, 0, true);
    printf("model loaded\n");
}

void terminate_stuff(){
    gsl_matrix_complex_free(w_Intergral_val);
    gsl_matrix_complex_free(w_Previous_error);

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