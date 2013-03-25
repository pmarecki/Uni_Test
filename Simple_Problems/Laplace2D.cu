#include <cstdio>
#include "stoper.h"

/**
 * Problem is ~N, ~bandwidth limited (large size). Async send/receive suggested.
 */

/**
 * Results at Quadro2000 (PM)
 */
//Computing Laplace A; A.width=128, Kernel.width=5
//H2D copy time=  0.60[msec]
//GPU exec time=  0.14[msec]
//CPUtime=  11.88[msec]
//Computing Laplace A; A.width=256, Kernel.width=5
//H2D copy time=  0.75[msec]
//GPU exec time=  0.25[msec]
//CPUtime=  49.52[msec]
//Computing Laplace A; A.width=512, Kernel.width=5
//H2D copy time=  1.40[msec]
//GPU exec time=  0.56[msec]
//CPUtime=  231.12[msec]
//Computing Laplace A; A.width=1024, Kernel.width=5
//H2D copy time=  8.26[msec]
//GPU exec time=  2.22[msec]
//CPUtime=  958.36[msec]
//Computing Laplace A; A.width=2048, Kernel.width=5
//H2D copy time=  8.78[msec]
//GPU exec time=  7.18[msec]
//CPUtime=  3853.60[msec]
//Computing Laplace A; A.width=4096, Kernel.width=5
//H2D copy time=  32.22[msec]
//GPU exec time=  28.02[msec]
//CPUtime=  15664.01[msec]

#define REP(i,n)  for(int i=0;i<(n);++i)
typedef uint32_t uint;


//Matrices being operated on
const int TSIZE = 32;               //all problems will be tiled
const int MSIZE = 8192;   //assume square matrices of this size
const int MSSIZE = MSIZE * MSIZE;

inline __host__ __device__ int pos(int x, int y){
  uint position = x + y*MSIZE;
  position %= MSIZE;
  return position;
}

const int kKernelSize = 5;
const int kBorder = kKernelSize / 2;    //separate algo for x<=border
inline __host__ __device__ int eK(int x, int y){
  uint position = x + y*kKernelSize;
  position %= kKernelSize;
  return position;
}

__global__ void Laplace(const float *A, const float *K, float *DDA);


int main(void) {
  // Allocation
  float *hA, *hDDA;              //allocated on host
  hA        = new float[MSSIZE];
  hDDA      = new float[MSSIZE];
  float LaplaceKernel[] = { 0,0,-1,0,0,  0,0,16,0,0,  -1,16,-60,16,-1,
                            0,0,-1,0,0,  0,0,16,0,0};
  REP(i, kKernelSize * kKernelSize)
    LaplaceKernel[2] /= 12.0;

  float *A, *Kernel, *DDA;                 //allocated on GPU
  cudaMalloc(&A, MSSIZE * 4);
  cudaMalloc(&Kernel, kKernelSize * kKernelSize * 4);
  cudaMalloc(&DDA, MSSIZE * 4);
  dim3 blocks(MSIZE / TSIZE, MSIZE / TSIZE);
  dim3 threads(TSIZE, TSIZE);   //1024 threads; standard configuration

  printf("Computing Laplace A; A.width=%i, Kernel.width=%i\n",
      MSIZE, kKernelSize);
  // Filling
  srand(12);
  for(int i=0; i<MSSIZE; ++i) {
    hA[i] = (rand() % 100 - 50) / 50.F;   //small mixed-size numbers
  }
  PosixStoper xx;
  cudaMemcpy(A, hA, MSSIZE * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(Kernel, LaplaceKernel, kKernelSize * kKernelSize * 4,
      cudaMemcpyHostToDevice);
  xx.Snap(); printf("H2D copy time=\t%3.2f[msec]\n",xx.LastDt()/1000);
  Laplace<<<blocks, threads>>>(A, Kernel, DDA);
  cudaDeviceSynchronize();
  xx.Snap(); printf("GPU exec time=\t%3.2f[msec]\n",xx.LastDt()/1000);
  cudaMemcpy(hDDA, DDA, MSSIZE * 4, cudaMemcpyDeviceToHost);

  float *htest = new float[MSSIZE];          //test results
  bzero(htest, MSSIZE * 4);
  PosixStoper yy;
  for(int x=0; x<MSIZE; ++x)
    for(int y=0; y<MSIZE; ++y)
    {
      for(int dx=-kBorder; dx<kBorder; ++dx)
        for(int dy=-kBorder; dy<kBorder; ++dy)
          htest[pos(x,y)] += hA[pos(x+dx,y+dy)] * LaplaceKernel[eK(dx,dy)];
    }
  yy.Snap(); printf("CPUtime=\t%3.2f[msec]\n",yy.LastDt()/1000);

  //Comparison
  int z=0;
  for(int x=0; x<MSIZE; ++x)
    for(int y=0; y<MSIZE; ++y)
      if (fabs(htest[pos(x,y)] - hDDA[pos(x,y)]) > 0.1) {
        ++z;
//        printf("ERR T[%i,%i]=%i\t C=%i\n", x,y, htest[x + y * kMatrixWidth],
//            hC[x + y * kMatrixWidth]);
    }
  printf("Err:%i\n",z);
  delete[] htest;
  cudaFree(A); cudaFree(Kernel); cudaFree(DDA);
  delete[] hA; delete hDDA;
}



/**
 * Each block computes a tile (bx=column,by=row) of C.
 * It must loop over a few tiles of A and B, and sum results.
 */
__global__ void Laplace(const float *A, const float *Kernel, float *DD) {
  // Tiles held in matrices sA, sB (SMEM), loaded by threads first.
  int bx = blockIdx.x;     //block-column in C  (column in B)
  int by = blockIdx.y;     //block-row    in C  (row    in A)
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int x = bx * blockDim.x + tx;
  int y = by * blockDim.y + ty;
  DD[pos(x,y)]=0;
  for(int dx=-kBorder; dx<kBorder; ++dx)
    for(int dy=-kBorder; dy<kBorder; ++dy)
      DD[pos(x,y)] += A[pos(x+dx,y+dy)] * Kernel[eK(dx,dy)];

}
