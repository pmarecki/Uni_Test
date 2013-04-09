#include "../STL_MACRO.h"

#define MMSIZE  MSIZE * MSIZE



__device__ inline int at(int x, int y) {
  x+= MSIZE;
  y+= MSIZE;
  x %= MSIZE;
  y %= MSIZE;
  int res = (x + y * MSIZE);
  //if (res<0 || res>=MMSIZE) printf("err");
  return res;
}

/**
 "psi"      : old wavefunction
 "n_psi"    : new wavefunction being constructed (at t=t+tau)
 "V"        : potential
 "norm"     : single place for each block where norm of the tile is written
*/


/**
 * Initialize potential "V", wavefunctions "psi, n_psi", and norm "norm".
 * Run on LARGE blocks: TSIZE x TSIZE.
 */
__global__ void Init(float *Psi, float *n_Psi, float *V, float* norm) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  Psi[at(x,y)] = 1.;
  n_Psi[at(x,y)] = 0;
  V[at(x,y)] = 0;

  if (x==0 || y==0 || x==MSIZE-1 || y==MSIZE-1)   //~potential well
    V[at(x,y)] = 1.;

  if (x < TILE_MX &&  y < TILE_MX)
      norm[x + y * TILE_MX]=0;
}

//  Run on 1 thread.
__global__ void Check(float *Psi, float *n_Psi) {
  float psisum = 0, npsisum = 0;
  REP(i,MSIZE)
    REP(j,MSIZE) {
      psisum += Psi[at(i,j)];
      npsisum += n_Psi[at(i,j)];
    }
  printf("Check: psi_sum=%3.2f(=%3.2f)\t npsi_sum=%3.2f(=0)\n", psisum,
      1.*MSIZE * MSIZE, npsisum);
}


// Run on blocks of size TSIZE x TSIZE
// Main Kernel
// Hamiltonian on psi Time:   64.107 [msec]  (TSIZE=6134)
// "Psi" is read-only; "n_Psi" is write-only.
__global__ void NPsi_HPsi(const float *Psi, float *n_Psi,
    const float *V, float* norm) {
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  /////////////
  // PART I
  // n_psi = Hamiltonian * psi
  // --->  -Laplace psi (periodic); Laplace = {{0,1,0}{1,-4,1}{0,1,0}}
  float psi = Psi[at(x,y)];       //local value of psi (at x,y)
  float n_psi = 0.;               //local value of n_psi (at x,y)
  n_psi = psi * 4;
  n_psi -= Psi[at(x,y-1)] + Psi[at(x,y+1)] + Psi[at(x-1,y)] + Psi[at(x+1,y)];
  n_psi /= hh;
  // ---> += V psi + g * psi^3
  n_psi += (V[at(x,y)] + g * psi * psi) * psi;
  // ---> n_psi = psi - tau * (H psi)
  n_psi *= -tau;
  n_psi += psi;
  n_Psi[at(x,y)] = n_psi;                 //store GMEM
   //////////////
   // PART II:
   // Computing norm-contribution from the whole block=tile.
   // Standard reduction algorithm.
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
      tile[ty][tx] += tile[ty][tx + size];
    size /= 2;
  }
  __syncthreads();
  //Now sum all tmp[i][0], for i=[0,TSIZE), and write to external "norm[]"
  if (tx==0 && ty==0)
  {
    float nrm = 0;
    REP(i, TSIZE)
      nrm += tile[i][0];
    norm[blockIdx.x + blockIdx.y * TILE_MX] = nrm;      //store GMEM
  }
}


//////////////////////////////////////////////////////
// Sums up numbers from TILE_MX * TILE_MX blocks (~256 * 256).
// Full norm stored into "norm[0]"
//////////////////////////////////////////////////////
// Run on <<<1,TILE_MX>>>; each thread summing up a single column of "norm[]".
// Used as finishing elements for kernels working on 2D objects,
// and summing things up (~normalization, or comutation of expectation values).
//////////////////////////////////////////////////////
__global__ void Reduce_finish(float* norm, int verbose) {
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
    REP(c, TILE_MX) {
      if (verbose)
        printf("{%3.4f}", norm[c]);
      column_norm += norm[c];
    }
    norm[0] = column_norm;
//    printf("\nFull norm sum:\t%3.4f\n", norm[0]);
    if (verbose)
    printf("\n");
  }

}


//Display content of (small) matrix "norm". Run on 1 thread.
__global__ void PrintNormMatrix(float* norm) {
  REP(i, TILE_MX) {
    REP(j, TILE_MX)
      printf("%3.3f\t", norm[i*TILE_MX+j]);
    printf("\n");
  }
}



// Copy normalized n_psi into psi;
// Run on blocks of size TSIZE x TSIZE
__global__ void NormalizePsi(float *psi, float *n_psi, float* norm) {
  float fact = sqrt((float)MMSIZE / norm[0]);
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  if (x==0 && y==0)
    printf("norm factor = %3.4f\n", fact);
  psi[at(x,y)] = n_psi[at(x,y)] * fact;
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
  n_psi = psi * 4;
  n_psi -= Psi[at(x,y-1)] + Psi[at(x,y+1)] + Psi[at(x-1,y)] + Psi[at(x+1,y)];
  n_psi /= hh;
  n_psi += (V[at(x,y)] + g * psi * psi) * psi;
  //<psi, H psi> = <psi, n_psi>
  n_psi *= psi;      //write into GMEM

//  ///////////
//  // Standard reduction algorithm.
//  // Summing up all contibutions from this tile
//
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
    norm[blockIdx.x + blockIdx.y * TILE_MX] = nrm;    //store GMEM
  }
  // Later: sum up the "norm[]" matrix by calling "Reduce_finish()"
}


__global__ void PrintEnergy(const float* norm) {
  printf("Energy = %3.4f\n", norm[0]);
}

__global__ void PrintPsi(const float* psi) {
  REP(y,10) {
    REP(i,10)
        printf("%3.1f ", psi[at(i,y)]);
    printf("\n");
  }
}






