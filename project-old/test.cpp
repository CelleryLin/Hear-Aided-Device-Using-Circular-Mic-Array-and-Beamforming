/*
  g++ -lpulse -lpulse-simple test.cpp
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

#include "AudioFile/AudioFile.h"

#define SAMPLE_RATE 44100
#define CH 6
using namespace std;

int main(){
    float sec=1;
    size_t bytesToWrite = SAMPLE_RATE*sec*CH;
    long int samples[bytesToWrite]={0};
    int error;

    pa_simple *rec_unit;
    pa_sample_spec recss;

    recss.format = PA_SAMPLE_S16NE;
    recss.rate = SAMPLE_RATE;
    recss.channels = CH;

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

        freopen("output.txt","w",stdout);
        pa_simple_read(rec_unit, samples, sizeof(samples), &error);
        //for(int i=0;i<bytesToWrite;i++) cout<<samples[i]<<"\n";


        //write to wav file
        int CHO=2;
        AudioFile <float> a;
        a.setNumChannels(CHO);
        a.setNumSamplesPerChannel(SAMPLE_RATE);
        a.setBitDepth(32);
        int k=0;
        for(int i=0;i<CHO;i++){
            for(int j=0;j<SAMPLE_RATE*sec;j++){
                a.samples[i][j]=(float)samples[j];
                k++;
                cout<<a.samples[i][j]<<"\n";
            }
        }

        a.save ("Output.wav", AudioFileFormat::Wave);
    return 0;
}