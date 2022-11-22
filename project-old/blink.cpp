#include <wiringPi.h>
#include "/home/pi/libapa102/lib/apa102.h"
#include <iostream>

int main(){
    //Initialize strip
    struct APA102* strip = APA102_Init(60);

    //Run animation
    struct APA102_Animation* anim = APA102_BlinkAnimation(strip, APA102_CreateFrame(31, 0xFF, 0x0, 0x0), 200);

    //Delay and kill
    delay(2000);
    APA102_KillAnimation(anim);
}