/*
  g++ -o -lpulse -lpulse-simple rec.out soundrec.cpp
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

#define SAMPLE_RATE 44100
#define CH 6
using namespace std;

int main(){
    float sec=1;
    size_t bytesToWrite = SAMPLE_RATE*sec*CH;
    //size_t bytesToWrite = 256*CH;
    long int samples[bytesToWrite]={0};
    int error;

    pa_simple *rec_unit,*play_unit;
    pa_sample_spec recss,plss;

    recss.format = PA_SAMPLE_S16NE;
    recss.rate = SAMPLE_RATE;
    recss.channels = CH;

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
                    NULL,
                    //"alsa_output.platform-bcm2835_audio.analog-stereo",               // Use the default device.
                    "Music",            // Description of our stream.
                    &plss,                // Our sample format.
                    NULL,               // Use default channel map
                    NULL,               // Use default buffering attributes.
                    NULL                // Ignore error code.
                    );

    freopen("output.txt","w",stdout);
    while(1){
        //cout<<"Listening...\n";
        pa_simple_read(rec_unit, samples, sizeof(samples), &error);
        //cout<<"Playing...\n";
        int k=0;
        long int newsamples[CH][(int)(SAMPLE_RATE*sec)]={0};
        for(int i=0;i<CH;i++){
            for(int j=0;j<SAMPLE_RATE*sec;j++){
                newsamples[i][j]=samples[k];
                k++;
                //cout<<newsamples[i][j]<<"\n";
            }
        }
        pa_simple_write(play_unit, newsamples[0], SAMPLE_RATE*sec*16, &error);
    }    
    return 0;
}