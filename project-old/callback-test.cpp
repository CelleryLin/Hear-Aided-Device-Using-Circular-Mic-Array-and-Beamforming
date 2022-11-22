/*
  g++ -lpulse -lpulse-simple -pthread callback-test.cpp
*/

#include<iostream>
#include<pulse/pulseaudio.h>
#include<pulse/simple.h>
#include<stdio.h>
#include<cmath>
#include<unistd.h>
#include<string.h>
#include<errno.h>
#include<thread>

#define SAMPLE_RATE 44100
#define CH 6

#define buffer_size 4096
#define bytesToWrite buffer_size*CH


using namespace std;


typedef int datatype;

void sinewave(float amp);
void playback(pa_simple *s, int *error);
void record(pa_simple *s, int *error);
void printval();
float freq=1000;
float sec=1;
//int buffer_size=4410;
//const size_t bytesToWrite = buffer_size*CH;
datatype samples[bytesToWrite]={0};

int main(){
  
  int error;
  
  pa_simple *play_unit, *rec_unit;
  pa_sample_spec plss, recss;
  recss.format = PA_SAMPLE_S32NE;
  recss.rate = SAMPLE_RATE;
  recss.channels = CH;
  
  plss.format = PA_SAMPLE_S32NE;
  plss.rate = SAMPLE_RATE;
  plss.channels = CH;

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

  //sinewave(1);
  /* ... and play it */
  //cout<<sizeof(samples)<<"\n";
  //cout<<sizeof(datatype)<<"\n";
  //cout<<bytesToWrite*sizeof(datatype)<<endl;
  thread rec_thread(record, rec_unit, &error);
  //thread printval_thread(printval);
  thread play_thread(playback, play_unit, &error);
  play_thread.join();
  rec_thread.join();
  //printval_thread.join();
  //playback(s, samples,&error);
  //pa_simple_write(s, samples, bytesToWrite, &error);
  pa_simple_free(play_unit);
  pa_simple_free(rec_unit);
  return 0;
}




void sinewave(float amp){
  //freopen("output.txt","w",stdout);
  //int j=0;
  for(int j=0;j<CH;j++){
    for(int i=0; i<buffer_size; ++i) {
      samples[i+j*buffer_size] = amp*pow(2,(sizeof(datatype)*8-1))*sin((2*M_PI*freq)*i/(SAMPLE_RATE*CH));
      //cout<<i<<": "<<samples[i]<<endl;
    }
  }
  //cout<<"hello!";
}

void playback(pa_simple *s, int *error){
  while(1){
    //datatype output_buffer[bytesToWrite]={0};
    //for(int i=0;i<bytesToWrite;i++){
    //  output_buffer[i]=samples[i];
    //}
    pa_simple_write(s, samples, sizeof(samples), error);
  }
}

void record(pa_simple *s, int *error){
  while(1){

    pa_simple_read(s, samples, bytesToWrite*sizeof(datatype), error);
  }
}

void printval(){
  int len=sizeof(samples) / sizeof(*samples);
  int micnum=1;
  float sum=0;
  while(1){
    sum=0;
    for(int i=0+micnum;i<len+micnum;i++){
      sum+=abs(samples[i]);
    }
    sum/=len;
    cout<<log(sum)<<endl;
  }
  
}