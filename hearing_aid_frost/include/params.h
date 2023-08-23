#define FRAME_BLOCK_LEN 1024
#define SAMPLING_RATE 48000
#define INPUT_CHANNEL 8
#define USED_CH 6
#define OUTPUT_CHANNEL 1
#define FILTER_LEN 6

# define M_PI 3.14159265358979323846
#define Radii 0.0463
#define SPEED 340

#define DOA 0 //rad
#define LEARNING_RATE 0.05

#define DELAY_STAGE_tau(n,doa) ((Radii/SPEED)*(cos(doa)-cos(doa-2*M_PI*n/USED_CH)))




void FIR_filter(float *input, float *output, float *w);
void signal_shift(float *input, float *output, int tau);
void set_u(float *fir);