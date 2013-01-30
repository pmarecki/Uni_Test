/**
 *  Test of Float2 type. 
 */

#include <iostream>
using namespace std;

#ifndef _CUDA_
 #define __device__
 #define __host__
#endif

#include "dsmath_def.h"     // definition of `float2` numbers
#include "dsmath_pi.h"	    // basic operators of `float2` numbers

int main(void) {
  float2 c;
  set(c, 12.);   //explicit way of setting a float2

  c = (float2) 12.;       //explicit constructor from double
  float2 d = (float2) 3;

  double dc = 12;
  double dd = 3;

  c = c / d;           // simple math in float2
  dc = dc / dd * c;    // conversion float2-->double (for `dd * c`)

  cout <<  (double)c << endl;
  cout <<  dc << endl;
}
