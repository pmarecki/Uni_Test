#include "../STL_MACRO.h"
#include "stoper.h"

/**
 * Simple Gross-Pitaevskii ground-state determination.
 *
 * Laplace: Abramowitz & Stegun 25.3.31 (p 885)
 */

//Matrices being operated on:
const int TSIZE    = 32;     //all problems will be tiled
const int MSIZE    = 256;   //main square matrices of this size; n * TSIZE
const int TILE_MX  = MSIZE / TSIZE;   // number of tiles

const float hh  = 0.0001;      //spatial const
const float g   = 0.01;        //nonlinearity constant
const float tau = 0.001;        //time-step

/**
 * ToDo:
 *  - simpler programs -- see what happens to uncertainties
 *  - trace all possible sources of race-conditions
 *  - find why results do not repeat
 */


#include "GP_Ground_Kernels.cu"



int main(void) {
  float *psi, *n_psi, *V, *norm;
  int size_f = sizeof(float);

  cudaMalloc(&psi, MSIZE * MSIZE * size_f);
  cudaMalloc(&n_psi, MSIZE * MSIZE * size_f);
  cudaMalloc(&V, MSIZE * MSIZE * size_f);
  cudaMalloc(&norm, TILE_MX * TILE_MX * size_f);

  dim3 blocks(TILE_MX, TILE_MX);
  dim3 threads(TSIZE, TSIZE);   //1024 threads; standard configuration

  printf("Gross-Pitaevskii ground state.\n");
  printf("Problem size (psi, n_psi, V): %iMB\n",
      MSIZE * MSIZE * size_f * 3 / 1024 / 1024);

  Init<<<blocks, threads>>>(psi, n_psi, V, norm);
  cudaDeviceSynchronize();


  for(int iter=0; iter<10; ++iter) {
    NPsi_HPsi<<<blocks, threads>>>(psi, n_psi, V, norm);
    cudaDeviceSynchronize();
    Reduce_finish<<<1, TILE_MX>>>(norm);
    cudaDeviceSynchronize();
    NormalizePsi<<<blocks, threads>>>(psi, n_psi, norm);
    cudaDeviceSynchronize();
    Energy<<<blocks, threads>>>(psi, V, norm);
    cudaDeviceSynchronize();
    Reduce_finish<<<1, TILE_MX>>>(norm);
    cudaDeviceSynchronize();
    PrintEnergy<<<1,1>>>(norm);
    cudaDeviceSynchronize();
  }


  cudaFree(psi);
  cudaFree(n_psi);
  cudaFree(V);
  cudaFree(norm);
}






