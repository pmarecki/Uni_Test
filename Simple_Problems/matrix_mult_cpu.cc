/**
 *  Test of tiled matrix multiplication on CPU. 
 *
 */ 

#include <cstdio>
#include <cstdlib>
#include "stoper.h"

const int TSIZE = 16;       //matrices divided into 32 x 32 tiles; CPU & GPU
const int MSIZE = 256;  //assume square matrices of this size
const int MSSIZE = MSIZE * MSIZE;

/**
 *  Scales quadratically only around kMatrixWidth ~128
 *  Later drops dramatically.
 */


#define REP(i,n) for(int i=0;i<(n);++i)


//implement tiling
void SimpleMultiply(int *hA, int *hB, int *htest) {
  for(int x=0; x<MSIZE; ++x)
    for(int y=0; y<MSIZE; ++y)
      for(int i=0; i<MSIZE; ++i) {
        htest[x + y *MSIZE] += hA[i + y * MSIZE] *
                                      hB[x + i * MSIZE];
      }
}

void LocalizedMultiply(int *hA, int *hB, int *htest) {
  int max_block = MSIZE / TSIZE;
  REP(x, max_block)       //x-tile in output
  REP(y, max_block)       //y-tile in output
  {
    int offx = x * TSIZE;
    int offy = y * TSIZE;
    REP(z, max_block)     //tile looped over
    {
      int offz = z * TSIZE;
      REP(tx, TSIZE)     //x in tile
      REP(ty, TSIZE)     //y in tile
      REP(tz, TSIZE) {   //z in tile        //introduce macros for [][]
          htest[tx + offx + (offy + ty) * MSIZE]
          += hA[tz + offz + (offy + ty) * MSIZE] *
             hB[tx + offx + (offz + tz) * MSIZE];
      }
    }
  }
}




int main(void) {
  // Allocation
  int *hA, *hB, *hC;  //allocated CPU
  hA = new int[MSSIZE];
  hB = new int[MSSIZE];
  hC = new int[MSSIZE];
  int *A, *B, *C;     //allocated on GPU
  printf("Matrix multiplication; width=height=%i\n",MSIZE);
  // Filling
  srand(12);
  for(int i=0; i<MSSIZE; ++i) {
    hA[i] =1;// rand() % 100 - 50;   //small mixed-size numbers
    hB[i] =1;// rand() % 100 - 50;   //small mixed-size numbers
  }
  // Test result on CPU
  int *htest = new int[MSSIZE];  //tested on CPU
  bzero(htest, MSSIZE * 4);  //clean it up
  PosixStoper yy;
//  SimpleMultiply(hA, hB, htest);
  LocalizedMultiply(hA, hB, htest);

  yy.Snap(); printf("CPUtime: %5.2f[msec]\n", yy.LastDt() / 1000);

  //Comparison
  int z = 0;
  for(int x=0; x<MSIZE; ++x)
    for(int y=0; y<MSIZE; ++y)
      if (htest[x + y * MSIZE] !=  MSIZE)
        ++z;
        //        printf("ERR T[%i,%i]=%i\n", x,y, htest[x + y * kMatrixWidth]);
  printf("Errs: %i\n",z);
  delete[] htest;
  delete[] hA; delete[] hB; delete hC;
}






