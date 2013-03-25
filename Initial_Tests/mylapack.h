/*
 * mylapack.h
 *
 *  Created on: Feb 19, 2013
 *      Author: fmarecki
 */

#ifndef MYLAPACK_H_
#define MYLAPACK_H_

extern "C" double dnrm2_(int *n, double *x, int* incx);  //Euclidean norm

//Solving linear system   "A * x = b", found "x" stored in "b".
extern "C" void dgesv_(int *n, int *nrhs, double *A, int *lda, int *ipiv,
    double *b, int *ldb, int *info);

//Copies x-->y
extern "C" void dcopy_(int *n, double *x, int *incx,
    double *y, int *incy);

// C = alpha * op(A) * OP(B) + beta * C
// transa=N|T|C  (eg. pointers to charN = 'N')
// So it works
extern "C" void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
    double *alpha, double *A, int *lda, double *B, int *ldb, double *beta,
    double *C, int *ldc);


//Computes eigenvalues
extern void dgeev_(char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info);



#endif /* MYLAPACK_H_ */
