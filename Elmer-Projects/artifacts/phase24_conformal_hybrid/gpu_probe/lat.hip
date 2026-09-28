#include <hip/hip_runtime.h>
#include <cstdio>
#include <chrono>
__global__ void k(double *x){ if(threadIdx.x==0 && x) x[0]+=1.0; }
int main(){
  double *d; hipMalloc(&d, 8*1024*1024*8);
  for(int i=0;i<100;i++) k<<<1,64>>>(d); hipDeviceSynchronize();
  auto t0=std::chrono::high_resolution_clock::now();
  for(int i=0;i<2000;i++) k<<<1,64>>>(d);
  hipDeviceSynchronize();
  auto t1=std::chrono::high_resolution_clock::now();
  for(int i=0;i<500;i++){ k<<<1,64>>>(d); hipDeviceSynchronize(); }
  auto t2=std::chrono::high_resolution_clock::now();
  double h; for(int i=0;i<500;i++) hipMemcpy(&h,d,8,hipMemcpyDeviceToHost);
  auto t3=std::chrono::high_resolution_clock::now();
  double *big; hipMalloc(&big, 64*1024*1024); hipMemset(big,0,64*1024*1024); hipDeviceSynchronize();
  auto t4=std::chrono::high_resolution_clock::now();
  for(int i=0;i<20;i++) hipMemcpy(big, d, 64*1024*1024 < 8*1024*1024*8 ? 64*1024*1024 : 8*1024*1024*8, hipMemcpyDeviceToDevice);
  hipDeviceSynchronize();
  auto t5=std::chrono::high_resolution_clock::now();
  auto us=[](auto a, auto b){ return std::chrono::duration<double,std::micro>(b-a).count(); };
  printf("async launch: %.1f us/launch\n", us(t0,t1)/2000);
  printf("launch+sync: %.1f us\n", us(t1,t2)/500);
  printf("8-byte D2H memcpy: %.1f us\n", us(t2,t3)/500);
  printf("D2D copy 64MB: %.1f GB/s\n", 20*64.0/1024/(us(t4,t5)*1e-6));
}
