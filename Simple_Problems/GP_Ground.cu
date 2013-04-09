#include "../STL_MACRO.h"
#include "../cudaErrors.h"
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

const float hh  = 0.001;      //spatial const
const float g   = 0.00;        //nonlinearity constant
const float tau = 0.01;        //time-step

const int verb = 1;


/**
 * ToDo:
 * - check results with analytic approach to 2d Schr"odinger problem
 * - wavefunction localized on edges -- wrong sign?
 * - print wavefunction 1d 2d; scripts for gnuplot
 */


#include "GP_Ground_Kernels.cu"



int main(void) {
  float *psi, *n_psi, *V, *norm;
  int size_f = sizeof(float);

  cudaMalloc(&psi,  MSIZE * MSIZE * size_f);
  cudaMalloc(&n_psi,MSIZE * MSIZE * size_f);
  cudaMalloc(&V,    MSIZE * MSIZE * size_f);
  cudaMalloc(&norm, TILE_MX * TILE_MX * size_f);

  //Standard configuration for kernels working on whole psi(x,y).
  dim3 blocks(TILE_MX, TILE_MX);
  dim3 threads(TSIZE, TSIZE);   //max 1024 threads if TSIZE=32


//  printf("Gross-Pitaevskii ground state.\n");
  printf("Problem size (psi, n_psi, V): %ikB\n",
      MSIZE * MSIZE * size_f * 3 / 1024 );
  printf("Large kernels launched with %iblocks of %ithreads\n",
      TILE_MX * TILE_MX, TSIZE*TSIZE);


  Init<<<blocks, threads>>>(psi, n_psi, V, norm);
  CUDA_CHK;
  cudaDeviceSynchronize();
//  Check<<<1,1>>>(psi, n_psi);


  for(int iter=0; iter<15; ++iter) {
    NPsi_HPsi<<<blocks, threads>>>(psi, n_psi, V, norm);
    PrintNormMatrix<<<1,1>>>(norm);
    Reduce_finish<<<1, TILE_MX>>>(norm, !verb);
    NormalizePsi<<<blocks, threads>>>(psi, n_psi, norm);
    cudaDeviceSynchronize();
    printf("psi=\n");
    PrintPsi<<<1,1>>>(psi);
    cudaDeviceSynchronize();
    printf("npsi=\n");
    PrintPsi<<<1,1>>>(n_psi);
    cudaDeviceSynchronize();
    Energy<<<blocks, threads>>>(psi, V, norm);
    Reduce_finish<<<1, TILE_MX>>>(norm, !verb);
    PrintEnergy<<<1,1>>>(norm);
    cudaDeviceSynchronize();
  }
  EEE();

  cudaFree(psi);
  cudaFree(n_psi);
  cudaFree(V);
  cudaFree(norm);
}






