#include <cstdio>
#include "stoper.h"

/**
 * This is slow. Should not be done by CPU>GPU>CPU. Just do it on CPU. 
 */

const int kTileSize = 16;       //matrices divided into 32 x 32 tiles; CPU & GPU
const int kMatrixWidth = 16;  //assume square matrices of this size
const int kSquareSize = kMatrixWidth * kMatrixWidth;

// Single kernel for (tiled, SMEM) matrix multiplication.
__global__ void Multiply(const int *A, const int *B, int *C);
// Straightforward matrix multiplication kernel. Doesn't use shared memory SMEM.
__global__ void SlowTranspose(const int *A, int *T);

__global__ void test() {
  printf("test\n");
}

int main(void) {
  // Allocation
  int *hA, *hT;  //allocated CPU
  hA = new int[kSquareSize];
  hT = new int[kSquareSize];
  int *A, *T;     //allocated on GPU
  cudaMalloc(&A, kSquareSize * 4); //size in Bytes
  cudaMalloc(&T, kSquareSize * 4); //size in Bytes
  printf("Matrix transpose; width=height=%i\n", kMatrixWidth);
  // Filling
  for(int i=0; i<kSquareSize; ++i)
    hA[i] = i;   //small mixed-size numbers
  // Copy to GPU
  PosixStoper xx;
  cudaMemcpy(A, hA, kSquareSize * 4, cudaMemcpyHostToDevice);
  // Transpose on GPU, and bring results back
  dim3 blocks(kMatrixWidth / kTileSize, kMatrixWidth / kTileSize);
  dim3 threads(kTileSize, kTileSize);
  //SlowTranspose<<<1, 1>>>(A, T);
  test<<<1,1>>>();
  cudaDeviceSynchronize();
  cudaMemcpy(hT, T, kSquareSize * 4, cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  xx.Snap(); printf("GPUtime=%3.4f[msec]\n",xx.LastDt()/1000);
  // Test result on CPU
  int *htest = new int[kSquareSize];  //tested on CPU
  bzero(htest, kSquareSize * 4);  //clean it up
  PosixStoper yy;
  int tmp=kMatrixWidth;
  for(int x=0; x<tmp; ++x)
    for(int y=0; y<tmp; ++y)
        htest[x + y *kMatrixWidth] = hA[y + x * kMatrixWidth];
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
  __shared__ float sA[kTileSize][kTileSize];       //"tile" matrices
  __shared__ float sB[kTileSize][kTileSize];
  // Loop over tiles, for each block in C seleceted by (bx,by)
  for(bk=0; bk < kMatrixWidth / kTileSize; ++bk) {
    //load matrices into SMEM
    sA[ty][tx] = A[(by * kTileSize + ty) * kMatrixWidth + (bk * kTileSize) +tx];
    sB[ty][tx] = B[(bk * kTileSize + ty) * kMatrixWidth + (bx * kTileSize) +tx];
    __syncthreads();
    // Multiple the tiles A * B --store--> C
    for(int k=0; k<kTileSize; ++k)
      Csub += sA[ty][k] * sB[k][tx];
    __syncthreads();
  }
  C[(by * kTileSize + ty) * kMatrixWidth + (bx * kTileSize + tx)] = Csub;
}


/**
 */
__global__ void SlowTranspose(const int *A, int *T) {
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  int posx = bx * kTileSize + tx;
  int posy = by * kTileSize + ty;
  T[posy * kMatrixWidth + posx] = A[posx * kMatrixWidth + posy];
  printf("%i",T[posy * kMatrixWidth + posx]);
}





