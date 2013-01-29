#include <stdio.h>
#include <cblas.h>  //necessary; compile with `-lblas`

//USE: https://developer.apple.com/library/mac/#documentation/Accelerate/
// Reference/BLAS_Ref/Reference/reference.html

extern double cblas_ddot(const int vec_size, const double *X,
    const int stride_inX, const double *Y, const int stride_inY);

extern void cblas_dcopy(const int N, const double *X, const int incX, double *Y,
    const int incY);


int main() {
  double a[3] = {1, 2, 3};
  double b[3] = {4, 5, 6};

  double dot_prod = cblas_ddot(3, a, 1, b, 1);
  printf(" The dot product is: %f \n", dot_prod);

  cblas_dcopy(3, a, 1, b, 1); //copies a->b
  for(int i=0; i<3; i++) printf("%3.0f,",b[i]);

  printf("\n");
  return 0;
};
