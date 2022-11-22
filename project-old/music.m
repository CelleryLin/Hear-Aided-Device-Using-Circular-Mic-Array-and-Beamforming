clc
clear all
pkg load communications;
format long %The data show that as long shaping scientific
doa=[20 60]/180*pi; %Direction of arrival
N=200;%Snapshots
w=[pi/4]';%Frequency
M=10;%Number of array elements
P=length(w); %The number of signal
lambda=150;%Wavelength
d=lambda/2;%Element spacing
snr=10;%SNA

[y, fs] = audioread("hilberttest.wav");
ds=0.0365;
dbit=fs*ds;

for ii=1:M
  x(ii,:)=hilbert(y(1+dbit*(ii-1):44100+dbit*(ii-1)));
endfor

x=x+awgn(x,snr);%Insert Gaussian white noise
R=x*x'; %Data covarivance matrix
[N,V]=eig(R); %Find the eigenvalues and eigenvectors of R
NN=N(:,1:M-P); %Estimate noise subspace


theta=-90:0.5:90; %Peak search
for ii=1:length(theta)
  SS=zeros(1,length(M));

  for jj=0:M-1
    SS(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
  end

  PP=SS*NN*NN'*SS';
  Pmusic(ii)=abs(1/ PP);
end

Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
plot(theta,Pmusic,'-k')
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title('DOA estimation based on MUSIC algorithm ')
grid on
