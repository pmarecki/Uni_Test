#include <cstdio>

__global__ void test() {
 printf("hello from TH: %i\n", threadIdx.x + blockIdx.x * blockDim.x);
}


int main(void) {
  test<<<2,2>>>();
  cudaDeviceSynchronize();
  return 0; 
}
