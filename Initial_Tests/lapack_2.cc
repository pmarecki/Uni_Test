#include <stddef.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cblas.h>
/**
 *  Compile via: "g++ test_lapack.cc -lblas -llapack"
 */

#include "mylapack.h"   //declarations of LAPACK functions.

int SIZE = 10;
int ONE = 1;


#define N 5
#define LDA N
#define LDVL N
#define LDVR N

int main(void) {
  int k=1;
  double A[SIZE * SIZE], b[SIZE * k], Acopy[SIZE * SIZE], bcopy[SIZE * SIZE];
  int ipiv[SIZE], info;
  int i;

  int ione = 1, size;
  char charN = 'N';
  double dpone = 1.e0, dmone = -1.e0;


  size = SIZE * k;
  dcopy_(&size, b, &ONE, bcopy, &ONE);

  size = SIZE * SIZE;
  dcopy_(&size, A, &ONE, Acopy, &ONE);
  dgesv_(&SIZE, &k, A, &SIZE, ipiv, b, &SIZE, &info);

  //matrix multiplication example
  dgemm_(&charN, &charN, &SIZE, &k, &SIZE, &dpone, Acopy,
      &SIZE, b, &SIZE, &dmone, bcopy, &SIZE);


  /* Locals */
  int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR,  lwork;
  double wkopt;
  double* work;
  /* Local arrays */
  double wr[N], wi[N], vl[LDVL * N], vr[LDVR * N];
  double a[LDA * N] = { -1.01, 3.98, 3.30, 4.43, 7.31, 0.86, 0.53, 8.26, 4.96,
      -6.43, -4.60, -7.04, -3.89, -7.66, -6.16, 3.31, 5.29, 8.20, -7.33, 2.47,
      -4.81, 3.55, -1.51, 6.18, 5.58 };
  /* Executable statements */
  printf(" DGEEV Example Program Results\n");
  /* Query and allocate the optimal workspace */
  lwork = -1;
  dgeev_("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &wkopt,
      &lwork, &info);
  lwork = (int) wkopt;
  work = (double*) malloc(lwork * sizeof(double));
  /* Solve eigenproblem */
  dgeev_("Vectors", "Vectors", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work,
      &lwork, &info);
  /* Check for convergence */
  if (info > 0) {
    printf("The algorithm failed to compute eigenvalues.\n");
    exit(1);
  }


  return 0;

}
