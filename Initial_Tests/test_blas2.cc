#include <cblas.h>  //necessary; compile with `-lblas`
#include "../STL_MACRO.h"

//USE: https://developer.apple.com/library/mac/#documentation/Accelerate/
// Reference/BLAS_Ref/Reference/reference.html

//extern double cblas_ddot(const int vec_size, const double *X,
//    const int stride_inX, const double *Y, const int stride_inY);
//
//extern void cblas_dcopy(const int N, const double *X, const int incX, double *Y,
//    const int incY);

void PrintDMA(double *ma, int sx, int sy) {
  REP(i,sx) {
    REP(j,sy)
      printf("%3.2f ", ma[i*sy+j]);
    EEE();
  }


}

int main() {
  double a[3] = {1, 2, 3};
  double b[3] = {4, 5, 6};

  double dot_prod = cblas_ddot(3, a, 1, b, 1);
  printf(" The dot product is: %f \n", dot_prod);

  cblas_dcopy(3, a, 1, b, 1); //copies a->b
  REP(i,3) printf("%3.0f,",b[i]); EEE();

  double sum = cblas_dasum(3, a, 1);
  printf("Sum of all elements is: %4.2f\n", sum);

  double *A = new double[9];
  double *B = A;
  double *C = new double[9];
  REP(i,9) A[i] = i % 3;
  PrintDMA(A,3,3);
  EEE();

  //C = alpha AB + beta C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
      3/*rowsAC*/, 3/*cols BC */, 3/*cols A rows B*/,
      1 /*alpha*/, A, 3 /*1st dim A*/, B, 3 /*1st dim B*/, 1.0 /*beta*/,
      C, 3 /*1st dim C*/);
  PrintDMA(C,3,3);

  delete[] A;
  delete[] C;
  printf("\n");
  return 0;
};
