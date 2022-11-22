#include <iostream>
#include <stdlib.h>
#include <thread>

using namespace std;

void thread1(long long int *res){
  
  *res=0;
  for(int i=1;i<=10000;i++){
    for(int j=1;j<=1000;j++){
      *res+=i*j;
      //cout<<"bbbbbbbbbbbb\n";
    }
  }
  
}

void thread2(long long int *res){
  *res=0;
  for(int i=10001;i<=20000;i++){
    for(int j=1;j<=1000;j++){
      *res+=i*j;
      //cout<<"hello\n";
    }
  }
}

void threadarr(int *arr){
  arr[0]=1;
}

int main(){
    long long int a,b,r;
    int arr[3000]={0};
    size_t bytesToWrite = 44100;
    int samples[bytesToWrite]={0};
    thread th1(thread1,&a);
    thread th2(thread2,&b);
    //thread th3(threadarr,samples);
    //threadarr(arr);
    th1.join();
    th2.join();
    //th3.join();
    //thread1(&a);
    //thread2(&b);
    r=a+b;
    cout<<r;
}