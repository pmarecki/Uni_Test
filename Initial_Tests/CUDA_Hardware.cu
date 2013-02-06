/**
 * @file  CUDA_Hardware.cu
 * @brief Report essential characteristics of GPU.
 *
 *
 */
#include <cstdio>
#include <cstdlib>
#include <iostream>
using namespace std;



int main(int argc, char** argv) {
  int count = 0;
  cudaGetDeviceCount(&count);
  printf("Report on GPU configuration (GPUs: %i).\n", count);
  for(int i=0; i<count; i++) {
    printf("Device:1/%i\n", count);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Global GPU Memory (GMEM) :\t%li(MB)\n", prop.totalGlobalMem / 1024 /1024);
    printf("Streaming Multiprocessors:\t%i\n", prop.multiProcessorCount);
    printf("L2 cache size            :\t%i(kB)\n", prop.l2CacheSize / 1024);
    int thsm = prop.maxThreadsPerMultiProcessor;
    int thbl = prop.maxThreadsPerBlock;
    printf("Max threads per block    :\t%i\n", thbl);
    if (thsm == 2048)
      printf("Compute capability:\t3.0 (Kepler)\n");
    else if (thsm == 1536)
      printf("Compute capability:\t2.0 (Fermi)\n");
    else
      printf("Compute capability:\t1.3\n");
  }
}


