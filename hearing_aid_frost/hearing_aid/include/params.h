#define FRAME_BLOCK_LEN 1200    // using 1200 due to fucking matlab
#define SAMPLING_RATE 48000     // using 48k due to fucking matlab
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

#define vad_threshold 20

#define snr_gate 1e-5          //treshold of snr
#define min_snr 1e-3
#define snr_start_freq 200.
#define snr_end_freq 1000.


#define GAIN 1
#define GAIN_BIG (GAIN*10)
#define EMA_RATE 0.01 // 0~1

#define wKp 0.7
#define wKi 0.0
#define wKd 0.

#define start_freq 600
#define end_freq 1000


gsl_complex double2complex(float re, float im);
float EMA(float a, float a0, float alpha);
float SNR(fftwf_complex a[],double m_freq, double bandwidth);
void PrintMat(gsl_matrix_complex *m, char info[]);
void autocorr(gsl_matrix_complex *Rxx, gsl_matrix_complex *m);
void invMat(gsl_matrix_complex *matrix,gsl_matrix_complex *inv);
gsl_complex eulers_formula(double x);
float Compress(float n);
void Compress_Mat(gsl_matrix_complex *matrix);
float energy_vad(fftwf_complex a[],float threshold);