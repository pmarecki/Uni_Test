#include <cstdio>
#include "stoper.h"

/**
 * This is slow. Should not be done by CPU>GPU>CPU. Just do it on CPU. 
 */

const int TSIZE = 16;       //matrices divided into 32 x 32 tiles; CPU & GPU
const int MSIZE = 16;  //assume square matrices of this size
const int MSSIZE = MSIZE * MSIZE;

// Single kernel for (tiled, SMEM) matrix multiplication.
__global__ void Multiply(const int *A, const int *B, int *C);
// Straightforward matrix multiplication kernel. Doesn't use shared memory SMEM.
__global__ void SlowTranspose(const int *A, int *T);

__global__ void funkcja() {
  printf("test\n");
}

int main(void) {
  // Allocation
  int *hA, *hT;  //allocated CPU
  hA = new int[MSSIZE];
  hT = new int[MSSIZE];
  int *A, *T;     //allocated on GPU
  cudaMalloc(&A, MSSIZE * 4); //size in Bytes
  cudaMalloc(&T, MSSIZE * 4); //size in Bytes
  printf("Matrix transpose; width=height=%i\n", MSIZE);
  // Filling
  for(int i=0; i<MSSIZE; ++i)
    hA[i] = i;   //small mixed-size numbers
  // Copy to GPU
  PosixStoper xx;
  cudaMemcpy(A, hA, MSSIZE * 4, cudaMemcpyHostToDevice);
  // Transpose on GPU, and bring results back
  dim3 blocks(MSIZE / TSIZE, MSIZE / TSIZE);
  dim3 threads(TSIZE, TSIZE);
  //SlowTranspose<<<1, 1>>>(A, T);
  funkcja<<<1,1>>>();
  cudaDeviceSynchronize();
  cudaMemcpy(hT, T, MSSIZE * 4, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  xx.Snap(); printf("GPUtime=%3.4f[msec]\n",xx.LastDt()/1000);
  // Test result on CPU
  int *htest = new int[MSSIZE];  //tested on CPU
  bzero(htest, MSSIZE * 4);  //clean it up
  PosixStoper yy;
  int tmp=MSIZE;
  for(int x=0; x<tmp; ++x)
    for(int y=0; y<tmp; ++y)
        htest[x + y *MSIZE] = hA[y + x * MSIZE];
  yy.Snap(); printf("CPUtime=%3.4f[msec]\n",yy.LastDt()/1000);

  //Comparison
  int z=0;
//  for(int x=0; x<kMatrixWidth; ++x)
//    for(int y=0; y<kMatrixWidth; ++y)
//      if (htest[x + y * kMatrixWidth] != hT[x + y * kMatrixWidth]) {
//        ++z;
//        printf("CPU:%i GPU:%i\n", htest[x + y * kMatrixWidth],
//            hT[x + y * kMatrixWidth]);
//      }

  printf("Err:%i\n",z);
  delete[] htest;
  cudaFree(A); cudaFree(T);
  delete[] hA; delete[] hT;
}





/**
 * Each block computes a tile (bx=column,by=row) of C.
 * It must loop over a few tiles of A and B, and sum results.
 */
__global__ void Multiply(const int *A, const int *B, int *C) {
  // Tiles held in matrices sA, sB (SMEM), loaded by threads first.
  int bx = blockIdx.x;     //block-column in C  (column in B)
  int by = blockIdx.y;     //block-row    in C  (row    in A)
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int bk;   //index for loop over block of _tiles_ in A (row) and B (column).
  int Csub = 0;  //Store locally data in loop; write to GMEM only once at end.
  __shared__ float sA[TSIZE][TSIZE];       //"tile" matrices
  __shared__ float sB[TSIZE][TSIZE];
  // Loop over tiles, for each block in C seleceted by (bx,by)
  for(bk=0; bk < MSIZE / TSIZE; ++bk) {
    //load matrices into SMEM
    sA[ty][tx] = A[(by * TSIZE + ty) * MSIZE + (bk * TSIZE) +tx];
    sB[ty][tx] = B[(bk * TSIZE + ty) * MSIZE + (bx * TSIZE) +tx];
    __syncthreads();
    // Multiple the tiles A * B --store--> C
    for(int k=0; k<TSIZE; ++k)
      Csub += sA[ty][k] * sB[k][tx];
    __syncthreads();
  }
  C[(by * TSIZE + ty) * MSIZE + (bx * TSIZE + tx)] = Csub;
}


/**
 */
__global__ void SlowTranspose(const int *A, int *T) {
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int posx = bx * TSIZE + tx;
  int posy = by * TSIZE + ty;
  T[posy * MSIZE + posx] = A[posx * MSIZE + posy];
  printf("%i",T[posy * MSIZE + posx]);
}





