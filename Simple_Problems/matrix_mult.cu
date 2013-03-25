#include <cstdio>
#include "stoper.h"

/**
 * Note: CPU code scales as N^3 only for sizes up to ~128 (linux, AMD,
 *   512kB cache L2). Later dramatically slower (almost as N^4).
 */

const int TSIZE = 32;       //matrices divided into 32 x 32 tiles; CPU & GPU
const int MSIZE = 256;  //assume square matrices of this size
const int MSSIZE = MSIZE * MSIZE;

// Single kernel for (tiled, SMEM) matrix multiplication.
__global__ void Multiply(const int *A, const int *B, int *C);
// Straightforward matrix multiplication kernel. Doesn't use shared memory SMEM.
__global__ void SlowMultiply(const int *A, const int *B, int *C);



int main(void) {
  // Allocation
  int *hA, *hB, *hC;  //allocated CPU
  hA = new int[MSSIZE];
  hB = new int[MSSIZE];
  hC = new int[MSSIZE];
  int *A, *B, *C;     //allocated on GPU
  cudaMalloc(&A, MSSIZE * 4); //size in Bytes
  cudaMalloc(&B, MSSIZE * 4); //size in Bytes
  cudaMalloc(&C, MSSIZE * 4); //size in Bytes
  printf("Matrix multiplication; width=height=%i\n",MSIZE);
  // Filling
  srand(12);
  for(int i=0; i<MSSIZE; ++i) {
    hA[i] = rand() % 100 - 50;   //small mixed-size numbers
    hB[i] = rand() % 100 - 50;   //small mixed-size numbers
  }
  // Copy to GPU
  cudaMemcpy(A, hA, MSSIZE * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(B, hB, MSSIZE * 4, cudaMemcpyHostToDevice);
  // Multiply on GPU, and bring results back
  dim3 blocks(MSIZE / TSIZE, MSIZE / TSIZE);
  dim3 threads(TSIZE, TSIZE);
  PosixStoper xx;
  Multiply<<<blocks, threads>>>(A, B, C);
  cudaMemcpy(hC, C, MSSIZE * 4, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  xx.Snap(); printf("GPUtime=%3.2f[msec]\n",xx.LastDt()/1000);
  // Test result on CPU
  int *htest = new int[MSSIZE];  //tested on CPU
  bzero(htest, MSSIZE * 4);  //clean it up
  PosixStoper yy;
  int tmp=MSIZE;
  for(int x=0; x<tmp; ++x)
    for(int y=0; y<tmp; ++y)
      for(int i=0; i<tmp; ++i)
        htest[x + y *MSIZE] += hA[i + y * MSIZE] *
                                      hB[x + i * MSIZE];
  yy.Snap(); printf("CPUtime=%3.2f[msec]\n",yy.LastDt()/1000);

  //Comparison
  int z=0;
  for(int x=0; x<MSIZE; ++x)
    for(int y=0; y<MSIZE; ++y)
      if (htest[x + y * MSIZE] != hC[x + y * MSIZE])
        ++z;
//        printf("ERR T[%i,%i]=%i\t C=%i\n", x,y, htest[x + y * kMatrixWidth],
//            hC[x + y * kMatrixWidth]);
  printf("Err:%i\n",z);
  delete[] htest;
  cudaFree(A); cudaFree(B); cudaFree(C);
  delete[] hA; delete[] hB; delete hC;
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
 * Each block computes a tile (bx=column, by=row) of C.
 * It must loop over a few tiles of A and B, and sum results.
 */
__global__ void SlowMultiply(const int *A, const int *B, int *C) {
  int bx = blockIdx.x;     //block-column in C  (column in B)
  int by = blockIdx.y;     //block-row    in C  (row    in A)
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int bk;   //index for loop over block of _tiles_ in A (row) and B (column).
  int Csub = 0;
  // Loop over tiles, for each block in C seleceted by (bx,by)
  for(bk=0; bk < MSIZE / TSIZE; ++bk) {
    for(int k=0; k<TSIZE; ++k) {  //loop over index in the tile
      Csub +=  A[(by * TSIZE + ty) * MSIZE + (bk * TSIZE) + k]
             * B[(bk * TSIZE + k) * MSIZE + (bx * TSIZE) + tx];
    }
    __syncthreads();
  }
  C[(by * TSIZE + ty) * MSIZE + (bx * TSIZE + tx)] = Csub;
}





