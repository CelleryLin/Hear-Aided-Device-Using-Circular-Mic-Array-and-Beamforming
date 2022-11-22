/*
  g++ -lpulse -lpulse-simple soundplay.cpp
*/

#include<iostream>
#include<pulse/pulseaudio.h>
#include<pulse/simple.h>
#include<stdio.h>
#include<cmath>

#define SAMPLE_RATE 44100
#define CH 2
using namespace std;

void sinewave(size_t bytesToWrite,int *samples,float amp);
void squarewave(size_t bytesToWrite,int *samples,float amp);
void playpack(pa_simple *s, int *data, size_t bytes, int *error);
float freq=400;
float sec=1;
int main(){
  size_t bytesToWrite = SAMPLE_RATE*sec*CH*4;
  int samples[bytesToWrite];
  int error;
  
  pa_simple *s;
  pa_sample_spec ss;
  
  ss.format = PA_SAMPLE_S32NE;
  ss.rate = SAMPLE_RATE;
  ss.channels = CH;

  s = pa_simple_new(NULL,               // Use the default server.
                    "Fooapp",           // Our application's name.
                    PA_STREAM_PLAYBACK,
                    NULL,
                    //"alsa_output.platform-bcm2835_audio.analog-stereo",               // Use the default device.
                    "Music",            // Description of our stream.
                    &ss,                // Our sample format.
                    NULL,               // Use default channel map
                    NULL,               // Use default buffering attributes.
                    NULL                // Ignore error code.
                    );






  sinewave(bytesToWrite,samples,1);
  /* ... and play it */
  cout<<sizeof(*samples);
  playpack(s, samples,bytesToWrite,&error);
  pa_simple_free(s);
  return 0;
}




void sinewave(size_t bytesToWrite,int *samples,float amp){
  for(int i=0; i<bytesToWrite; ++i) {
      samples[i] = amp*pow(2,(sizeof(*samples)*8-1))*sin((2*M_PI*freq)*i/(SAMPLE_RATE*CH));
  }
}

void playpack(pa_simple *s, int *samples, size_t bytesToWrite, int *error){
  while(1){
    pa_simple_write(s, samples, bytesToWrite, error);
  }
}