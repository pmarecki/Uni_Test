/**
 *  Test of Float2 type. 
 */

#include <cstdio>
#include <cmath>
using namespace std;

#define _XGPU_
#include "FF.h"     // definition of `float2` numbers

/**
 * TODO: _) write explicit rigorous tests
 *       _) convert everything back to __device__ __host__
 */

// Basic creation and conversions<--->double;
//void TEST1() {
//  int i;
//}



__global__ void aaa() {
  printf("hahass\n");
}

int main(void) {
  double x = sin(1.02);
  FF x2= (FF)x;   //explicit constructor from double
//  x2 = (Float2)x;         // "x2=x;" would not compile
  //  x = x2;                 // this is allowed (conversion to double)

  x2 = x2 * x2 * x2 * x2;

//  Float2 z;
//  z=x2;

  double y;// = (double)x2;
  set(y, x2);

  printf("(f2)%2.16f : (d)%2.16f; diff= %2.16f\n", y, x*x*x*x, y-x*x*x*x);


}
