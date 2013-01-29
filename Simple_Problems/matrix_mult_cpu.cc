/**
 *  Test of tiled matrix multiplication on CPU. 
 *
 */ 

#include <cstdio>
#include <cstdlib>
#include "stoper.h"

const int kTileSize = 16;       //matrices divided into 32 x 32 tiles; CPU & GPU
const int kMatrixWidth = 256;  //assume square matrices of this size
const int kSquareSize = kMatrixWidth * kMatrixWidth;

/**
 *  Scales quadratically only around kMatrixWidth ~128
 *  Later drops dramatically.
 */


#define REP(i,n) for(int i=0;i<(n);++i)


//implement tiling
void SimpleMultiply(int *hA, int *hB, int *htest) {
  for(int x=0; x<kMatrixWidth; ++x)
    for(int y=0; y<kMatrixWidth; ++y)
      for(int i=0; i<kMatrixWidth; ++i) {
        htest[x + y *kMatrixWidth] += hA[i + y * kMatrixWidth] *
                                      hB[x + i * kMatrixWidth];
      }
}

void LocalizedMultiply(int *hA, int *hB, int *htest) {
  int max_block = kMatrixWidth / kTileSize;
  REP(x, max_block)       //x-tile in output
  REP(y, max_block)       //y-tile in output
  {
    int offx = x * kTileSize;
    int offy = y * kTileSize;
    REP(z, max_block)     //tile looped over
    {
      int offz = z * kTileSize;
      REP(tx, kTileSize)     //x in tile
      REP(ty, kTileSize)     //y in tile
      REP(tz, kTileSize) {   //z in tile        //introduce macros for [][]
          htest[tx + offx + (offy + ty) * kMatrixWidth]
          += hA[tz + offz + (offy + ty) * kMatrixWidth] *
             hB[tx + offx + (offz + tz) * kMatrixWidth];
      }
    }
  }
}




int main(void) {
  // Allocation
  int *hA, *hB, *hC;  //allocated CPU
  hA = new int[kSquareSize];
  hB = new int[kSquareSize];
  hC = new int[kSquareSize];
  int *A, *B, *C;     //allocated on GPU
  printf("Matrix multiplication; width=height=%i\n",kMatrixWidth);
  // Filling
  srand(12);
  for(int i=0; i<kSquareSize; ++i) {
    hA[i] =1;// rand() % 100 - 50;   //small mixed-size numbers
    hB[i] =1;// rand() % 100 - 50;   //small mixed-size numbers
  }
  // Test result on CPU
  int *htest = new int[kSquareSize];  //tested on CPU
  bzero(htest, kSquareSize * 4);  //clean it up
  PosixStoper yy;
//  SimpleMultiply(hA, hB, htest);
  LocalizedMultiply(hA, hB, htest);

  yy.Snap(); printf("CPUtime: %5.2f[msec]\n", yy.LastDt() / 1000);

  //Comparison
  int z = 0;
  for(int x=0; x<kMatrixWidth; ++x)
    for(int y=0; y<kMatrixWidth; ++y)
      if (htest[x + y * kMatrixWidth] !=  kMatrixWidth)
        ++z;
        //        printf("ERR T[%i,%i]=%i\n", x,y, htest[x + y * kMatrixWidth]);
  printf("Errs: %i\n",z);
  delete[] htest;
  delete[] hA; delete[] hB; delete hC;
}






