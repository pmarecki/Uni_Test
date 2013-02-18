#include "STL_MACRO.h"
#include "nrutil_nr_cuda.h"

/**
 * "nrutil_nr_cuda.h" are modified headers of NR's, with accessors
 * to raw matrices. (Needed for movement of data to/from GPU.)
 *
 * Program tests this functionality.
 */



#define SIZE 10

template <typename T>
void Print2D(T* ptr, int sizex, int sizey) {
  REP(i, sizex) {
    REP(j, sizey)
      printf("%2i,", ptr[i*sizey+j]);
    printf("\b\n");
  }
}


template <typename T>
__device__
void DPrint2D(T* ptr, int sizex, int sizey) {
  REP(i, sizex) {
    REP(j, sizey)
      printf("%2i,", ptr[i*sizey+j]);
    printf("\b\n");
  }
}


__global__
void blabla(int* ptr, int sizex, int sizey) {
  if (threadIdx.x==0) {
    DPrint2D(ptr, sizex, sizey);
  }
}



int main() {
  NRMat3d<int> mat(SIZE, SIZE, 1);
  REP(i, SIZE)
    REP(j, SIZE)
    mat[i][j][0] = i + j;

  Print2D(mat.data(), mat.dim1(), mat.dim2());

  NRMat<int> cc(SIZE,SIZE);
  Print2D(cc.data(), cc.nrows(), cc.ncols());

  //now: move the matrix to GPU and print it from there
  int *p;
  cudaMalloc(&p,4*SIZE*SIZE);
  cudaMemcpy(p, mat.data(), 4*SIZE*SIZE, cudaMemcpyHostToDevice);
  blabla<<<1,1>>>(p,SIZE, SIZE);
  cudaFree(p);

}
