#include "STL_MACRO.h"
#include "NR_MAT.h"

const int kSize = 10;

template <typename T>
void PM(NRMat3d<T> &mat, const char* what) {
  printf("Matrix %s\n", what);
  REP(i, mat.dim1()) {
    REP(j, mat.dim2())
      printf("%2i,", mat[i][j][0]);
    printf("\b\n");
  }
}

/**
 * Further goals:
 *
 * Define simple Laplace operator /registered sizes/ applicable trivially to
 * the GrossPitaevskii code (appears 3 times there), in CPU and GPU case.
 *
 * See how this can be done for FF numbers.
 */



int main() {
  NRMat3d<int> mat(kSize, kSize, 1);
  REP(i, kSize) REP(j, kSize)
    mat[i][j][0] = i + j;

  PM(mat, "A");

  int *raw_matrix = mat.RawMatrix2D();
  int sizex = mat.dim1();
  int sizey = mat.dim2();

  REP(i, sizex*sizey)
   printf("%i, ", raw_matrix[i]);
  EEE();



}
