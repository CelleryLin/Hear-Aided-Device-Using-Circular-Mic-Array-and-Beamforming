/*
  g++ -lpulse -lpulse-simple -pthread callback-test.cpp
*/

#include<iostream>
#include<pulse/pulseaudio.h>
#include<pulse/simple.h>
#include<pulse/error.h>
#include<stdio.h>
#include<cmath>
#include<unistd.h>
#include<string.h>
#include<errno.h>
#include<thread>

#define SAMPLE_RATE 44100
#define CH 6
using namespace std;

void sinewave(size_t bytesToWrite,int *samples,float amp);
void playback(pa_simple *s, int *data, size_t bytes, int *error);
void record(pa_simple *s, int *data, size_t bytes, int *error);

int main(){
    float sec=1;
    const size_t Input_Buffer = 1024;
    int samples[Input_Buffer*CH]={0};
    int error;

    pa_simple *rec_unit,*play_unit;
    pa_sample_spec recss,plss;

    recss.format = PA_SAMPLE_S16NE;
    recss.rate = SAMPLE_RATE;
    recss.channels = 6;

    plss.format = PA_SAMPLE_S16NE;
    plss.rate = SAMPLE_RATE;
    plss.channels = 6;
    rec_unit = pa_simple_new(NULL,               // Use the default server.
                      "Fooapp",           // Our application's name.
                      PA_STREAM_RECORD,
                      //NULL,
                      "alsa_input.platform-soc_sound.multichannel-input",               // Use the default device.
                      "Music",            // Description of our stream.
                      &recss,                // Our sample format.
                      NULL,               // Use default channel map
                      NULL,               // Use default buffering attributes.
                      NULL                // Ignore error code.
                      );

    play_unit = pa_simple_new(NULL,               // Use the default server.
                      "Fooapp",           // Our application's name.
                      PA_STREAM_PLAYBACK,
                      //NULL,
                      "alsa_output.platform-bcm2835_audio.analog-stereo",               // Use the default device.
                      "Music",            // Description of our stream.
                      &plss,                // Our sample format.
                      NULL,               // Use default channel map
                      NULL,               // Use default buffering attributes.
                      NULL                // Ignore error code.
                      );

    //sinewave(Input_Buffer,samples,1);
    //thread rec_thread(record,rec_unit, samples, Input_Buffer*6, &error);
    sinewave(Input_Buffer,samples,1);
    thread play_thread(playback,play_unit, samples, Input_Buffer*6, &error);
    //playback(play_unit, samples, Input_Buffer, &error);
    //rec_thread.join();
    //play_thread.join();
    
    return 0;
}


void record(pa_simple *s, int *data, size_t bytes, int *error){
    cout<<"hello";
    while(true){
        pa_simple_read(s, data, bytes, error);
    }
}

void playback(pa_simple *s, int *data, size_t bytes, int *error){
    while(true){
        pa_simple_write(s, data, bytes, error);
    }
}

void sinewave(size_t bytesToWrite,int *samples,float amp){
  for(int j=0;j<CH;j++){
    for(int i=0; i<bytesToWrite; ++i) {
    samples[i+bytesToWrite*j] = amp*pow(2,(sizeof(*samples)*8-1))*sin((2*M_PI*440)*i/(SAMPLE_RATE*CH));
    }
  }
}