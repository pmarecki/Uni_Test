
#include "../STL_MACRO.h"


//Helper: position in full matrix (periodicity implemented)
inline __host__ __device__
int at(int x, int y){
  uint position = MSIZE + x + y * MSIZE;      //also OK for negative "x"
  position %= MSIZE;
  return position;
}

/**
 "psi"      : old wavefunction
 "n_psi"    : new wavefunction being constructed (at t=t+tau)
 "V"        : potential
 "norm"     : single place for each block where norm of the tile is written
*/


/**
 * Initialize potential "V", wavefunctions "psi, n_psi", and norm "norm".
 * Run on blocks of size TSIZE x TSIZE.
 */
__global__ void Init(float *psi, float *n_psi, float *V, float* norm) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  psi[at(x,y)] = 1;
  n_psi[at(x,y)] = 0;

  if (x==0 || y==0 || x==MSIZE-1 || y==MSIZE-1)   //~potential well
    V[at(x,y)] = 5;
  else
    V[at(x,y)] = 0;

  if (x < TILE_MX &&  y < TILE_MX)
      norm[x + y * TILE_MX]=0;
}






// Run on blocks of size TSIZE x TSIZE
// Main Kernel
// Hamiltonian on psi Time:   64.107 [msec]  (TSIZE=6134)
// "Psi" is read-only; "n_Psi" is write-only.
__global__ void NPsi_HPsi(const float *Psi, float *n_Psi, float *V, float* norm) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  /*
   * Part I:
   * n_psi = Hamiltonian * psi
   * (takes up half of the time of this kernel)
   */
  // --->  -Laplace psi (periodic)
  float psi = Psi[at(x,y)];       //local value of psi (at x,y)
  float n_psi = 0;                //local value of n_psi (at x,y)


  n_psi -= psi * 4;
  n_psi += Psi[at(x,y-1)]+ Psi[at(x,y+1)]+ Psi[at(x-1,y)]+ Psi[at(x+1,y)];
  n_psi /= hh;
  // ---> += V psi + g * psi^3
  n_psi += (V[at(x,y)] + g * psi * psi) * psi;

  ////   new_p = psi + tau * (H psi)
  n_psi *= tau;
  n_psi += psi;

  ////TEST
//  n_psi = 1;   //norm should give TSIZE * TSIZE

  n_Psi[at(x,y)] = n_psi;      //write into GMEM


  /*
   * Part II:
   * Computing norm-contribution from the whole block=tile.
   * Standard reduction algorithm.
   * (Takes up half of time of this kernel.)
   */
  __syncthreads();
  // Sum up the norm-contibutions from each element.
  __shared__ float tile[TSIZE][TSIZE];
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  tile[ty][tx] = n_psi * n_psi;    //load into SMEM

  // First sum up elements of rows tmp[i][j] --> into tmp[i][0]
  int size = TSIZE / 2;
  while(size > 0) {
    __syncthreads();
    if (tx < size)
      tile[ty][tx] += tile[ty][tx+size];
    size /= 2;
  }
  __syncthreads();
  //Now sum all tmp[i][0], for i=[0,TSIZE), and write to external "norm[]"
  if (tx==0 && ty==0) {
    float nrm = 0;
    REP(i, TSIZE)
      nrm += tile[i][0];
    //Write at single position in norm, determined by "blockIdx".
    norm[blockIdx.x + blockIdx.y * TILE_MX] = nrm;
//    printf("norm of block (%i,%i)=%3.2f\n", blockIdx.x, blockIdx.y, nrm);
  }

  //Later:
  //Reduce_finish()   //global norm
  //NormalizePsi()
}

//////////////////////////////////////////////////////
// Sums up numbers from TILE_MX * TILE_MX blocks (~256 * 256).
// Full norm stored into "norm[0]"
//////////////////////////////////////////////////////
// Run on <<<1,TILE_MX>>>; each thread summing up a single column of "norm[]"
// Used as finishing elements for kernels working on 2D objects,
// and summing things up (~normalization, or comutation of expectation values).
//////////////////////////////////////////////////////
__global__ void Reduce_finish(float* norm) {
  int x = threadIdx.x;   // in [0, TILE_MX)
  //Each thread sums up single column in matrix "norm[]"
  float column_norm = 0;
  REP(i, TILE_MX)
    column_norm += norm[x + i * TILE_MX];
  norm[x] = column_norm;
  __syncthreads();
  //All norms now in the single, first row.

  //Finally sum them all up.
  if (x==0) {
    column_norm = 0;
    REP(c, TILE_MX)
     column_norm += norm[c];
    norm[0] = column_norm;
    printf("Reduce --> %3.2f\n", norm[0]);
  }
}


// Copy normalized n_psi into psi;
// SMALL kernel. Run on blocks of size TSIZE x TSIZE
__global__ void NormalizePsi(float *psi, float *n_psi, float* norm) {
  float divisor = 1/sqrt(norm[0]);
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  psi[at(x,y)] = n_psi[at(x,y)] * divisor;
}


///////////////////////////////////////////////////
// Computation of <psi, H psi>
// Vector "norm" used as storage.
// Psi assumed normalized.
///////////////////////////////////////////////////

// Run on blocks of size TSIZE x TSIZE
__global__ void Energy(const float *Psi, const float *V, float* norm) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  // --->  -Laplace psi
  float psi = Psi[at(x,y)];       //local value of psi (at x,y)
  float n_psi = 0;                //local value of n_psi (at x,y)

  n_psi -= psi * 4;
  n_psi += Psi[at(x,y-1)] + Psi[at(x,y+1)] + Psi[at(x-1,y)] + Psi[at(x+1,y)];
  n_psi /= hh;
  // ---> += V psi + g * psi^3
  n_psi += (V[at(x,y)] + g * psi * psi) * psi;

  //<psi, H psi> = <psi, n_psi>
  n_psi *= psi;      //write into GMEM

  ///////////
  // Standard reduction algorithm.
  // Summing up all contibutions from this tile

  __syncthreads();
  // Sum up the norm-contibutions from each element.
  __shared__ float tile[TSIZE][TSIZE];
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  tile[ty][tx] = n_psi;       //load into SMEM

  // First sum up elements of rows tmp[i][j] --> into tmp[i][0]
  int size = TSIZE / 2;
  while(size > 0) {
    __syncthreads();
    if (tx < size)
      tile[ty][tx] += tile[ty][tx+size];
    size /= 2;
  }
  __syncthreads();
  //Now sum all tmp[i][0], for i=[0,TSIZE), and write to external "norm[]"
  if (tx==0 && ty==0) {
    float nrm = 0;
    REP(i, TSIZE)
      nrm += tile[i][0];
    //Write at single position in norm, determined by "blockIdx".
    norm[blockIdx.x + blockIdx.y * TILE_MX] = nrm;
    if (blockIdx.x == 0 && blockIdx.y == 0)
    printf("energy of block (%i,%i)=%3.2f\n", blockIdx.x, blockIdx.y, nrm);
  }
  // Later: sum up the "norm[]" matrix by calling "Reduce_finish()"
}

__global__ void PrintEnergy(const float* norm) {
  printf("Energy = %3.2f\n", norm[0]);
}






