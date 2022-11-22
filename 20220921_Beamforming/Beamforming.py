'''
20220921 Beamforming

Using MVDR Method to filter out the unwanted sound.
Using MUSIC Algorithm to estimate DOA dynamically (every 5 secs) 


#TODO
- Get rid of annoying beep sound (maybe rewrite the program using C++)
- Smooth the beamforming weight vector (w) changes (maybe use PID control)
- Using wide bandwidth beamforming

'''


import pyaudio
import time
import numpy as np
import scipy.linalg as LA
import struct
import threading as td
import multiprocessing as mp
import sys
from scipy.signal import hilbert
import scipy.signal as ss
from pyargus import directionEstimation as de
#import matplotlib.pyplot as plt

#np.set_printoptions(threshold=np.inf)

FORMAT = pyaudio.paInt32
CHANNELS = 8

ELEMENT=6
SOURCE=3
Radii=46.3e-3
LDA=343/2000; #wavelength (440hz for v=343m/s)
phi=np.radians(90)
Angles = np.arange(0,360,1)
Pmusic = np.zeros(len(Angles))
scanning_vectors=de.gen_uca_scanning_vectors(ELEMENT, Radii, Angles)

CHUNK = 2**10
RATE = 44100
RESPEAKER_INDEX = 1

input_buffer=np.zeros([CHANNELS,CHUNK])

p = pyaudio.PyAudio()

def MVDR(input_buffer,doa):
    b=input_buffer[0:ELEMENT]
    CovMat = de.corr_matrix_estimate(b.T, imp="fast")
    doa_processing=doa
    if (LA.det(CovMat)!=0 and len(doa_processing)!=0):        
        invR=LA.inv(CovMat)
        temp=np.sin(phi)*np.cos(np.radians(doa_processing[0])-2*np.arange(0,ELEMENT)*np.pi/ELEMENT)
        AA=np.exp(-1j*2*np.pi*Radii*temp/LDA)
        ww=np.matmul(np.matmul(AA,invR),AA.conj().transpose())
        #print("ww : {}".format(ww))
        w=np.matmul(invR,AA.conj().transpose())*(1/ww)
        
        return w
        #output_buffer=np.matmul(w.conj().transpose(),b)
        #print(output_buffer.shape)
        #return np.real(output_buffer)
        
        

def MU_Processing():
    global input_buffer, doa
    #print(input_buffer)
    b=hilbert(input_buffer[0:ELEMENT])
    CovMat = de.corr_matrix_estimate(b.T, imp="fast")
    MUSIC=de.DOA_MUSIC(CovMat, scanning_vectors, signal_dimension=SOURCE)
    #print(CovMat)
    Pmusic=np.abs(1/MUSIC)
    Pmusic=np.log10(Pmusic/np.amax(Pmusic))
    doa_index,_=ss.find_peaks(Pmusic,height=-0.8)
    doa=Angles[doa_index]
    
    #print(td.currentThread())
    return doa

def printToFile(a):     #Print something to ./filename.txt
    original_stdout = sys.stdout
    with open('./filename.txt', 'w') as f:
        sys.stdout = f
        print(a)
        sys.stdout = original_stdout

def makeArray(a):
    """
    transpose the original data.

    return data:
    [[mic1 t1, mic1 t2, mic1 t3 ... ],
    [mic1 t1, mic1 t2, mic1 t3 ... ],
                .
                .
                .
    [micN t1, micN t2, micN t3 ... ]]

    """
    b=np.zeros([CHANNELS,CHUNK])
    for i in range(CHANNELS):
        b[i][:]=a[i::CHANNELS]
    return b

def mix(a):
    b=np.zeros(CHUNK)
    c=[]
    d=[]
    for i in range(CHANNELS):
        b+=a[i::CHANNELS]
    b/=CHANNELS

    #b=a[0::CHANNELS]

    for i in range(CHANNELS):
        c.append(b)

    d=np.reshape(np.array(c).T,(CHANNELS*CHUNK))    
    return d

def monoToMu(a,ch,amp):
    a*=amp
    c=[]
    d=[]
    for i in range(ch):
        c.append(a)

    d=np.reshape(np.array(c).T,(CHANNELS*CHUNK))    
    return d


def callback(in_data, frame_count, time_info, status):
    """ 
    indata: 1D byte array received from mics.
    [mic1 t1, mic2 t1, mic3 t1 ... mic1 t2, mic2 t2, mic3 t2 ... micN tM]
    where N = CHANNELS, M=CHUNK.
    
    To make in_data calculable, convert it to nparray and transpose it.
    """
    #global CovMat
    global input_buffer, doa, w
    a=[]
    x=np.array(struct.unpack(str(frame_count*CHANNELS)+'i',in_data),dtype='int32')  #convert to nparray
    a=makeArray(x)
    input_buffer=a
    y=np.real((np.matmul(w.conj().transpose(),a[0:ELEMENT])))

    if (y.shape == (CHUNK,)):
        #print("YES")

        yy=monoToMu(y,CHANNELS,2)
        #yy=a
        #yy=mix(x)
        yy=yy.astype(dtype="int32")
        #print(y)
        out_data=yy.tobytes()
        return (out_data, pyaudio.paContinue)
    else:
        print("Mute")
        yy=np.zeros(CHUNK*CHANNELS)
        yy=yy.astype(dtype="int32")
        out_data=yy
        return (out_data, pyaudio.paContinue)

        #print("Bypassed")
        #yy=mix(x)
        #yy=yy.astype(dtype="int32")
        ##print(y)
        #out_data=yy.tobytes()
        #return (out_data, pyaudio.paContinue)


def main():
    global input_buffer, doa, w
    doa=[0]
    w=np.zeros([ELEMENT,ELEMENT])
    stream = p.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=RATE,
                    input=True,
                    output=True,
                    stream_callback=callback,
                    frames_per_buffer=CHUNK,
                    input_device_index=RESPEAKER_INDEX,
                    output_device_index=RESPEAKER_INDEX
                )


    
    #MUtd = td.Thread(target=MU_Processing,args=(1,),daemon=True)
    #MUtd.start()
    stream.start_stream()
    #BFtd.start()
    #MVDR(1)
    #heavycalctest(0)
    while stream.is_active():
        MU_Processing()
        w=MVDR(input_buffer,doa)
        print(doa)
        print(w)
        time.sleep(5)

    stream.stop_stream()
    stream.close()

    p.terminate()

if ( __name__== "__main__"):
    main()