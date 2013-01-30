#ifndef _FLOAT2_H_
#define _FLOAT2_H_

#include <inttypes.h>




//TODO be sure not to allow FMAD in multiplication
//CUDA might insert FMAD here -- find way out
#ifdef XGPU
__host__ __device__
#endif
float multiply(float a, float b)
{
  return a * b;
}

/**
 * Float2's are:
 *   explicitly constructed from `double`,
 *   have `operator=`,
 *   have `operator double`
 */
struct FF {
  float x, y;
//
  FF() : x(0), y(0) {}

//
  FF(double d) {
    x = (float) d;
    y = (float) (d - x);
  }
//
  FF& operator=(const FF &that) {
    x = that.x;
    y = that.y;
    return *this;
  }

  //will be augumented to "explicit operator double()", when c++11 arrives
  operator double() const {
    return (double) x + (double) y;
  }
};



/////////////////////////////
///  Arithmetic for Float2
///
// Double single functions based on DSFUN90 package:
// http://crd.lbl.gov/~dhbailey/mpdist/index.html
// Partially adapted from NVIDIA's CUDA Mandelbrot demo

// Set functions for agnostic arithmetic
 void set(FF &a, double b);
 void set(float  &a, double b);
 void set(double &a, FF b) {a = (double)b.x + b.y;}

  void set(FF &a, const float  b);
  void set(float  &a, const FF b);
  void set(FF &a, const FF b);
  void set(float  &a, const float  b);

// Arithmetic operators
 FF operator-(const FF a);
 FF operator+(const FF a, const FF b);
 FF operator-(const FF a, const FF b);
 FF operator*(const FF a, const FF b);
 FF operator/(const FF a, const FF b);

// Note that the `operator/` potentially makes use of the FMAD
// instruction and may lose some accuracy as a result.

// This function sets the DS number A equal to the double
// precision floating point number B.
 inline void set(FF &a, double b) {
  a.x = (float)b;
  a.y = (float)(b - a.x);
}

// Store a (truncated) double in a float
 inline void set(float& a, double b) {
  a = (float)b;
}

// Store a float into a double single
  inline void set(FF &a, const float b) {
  a.x = b;
  a.y = 0;
}

// Store the hi word of the double single in a float
  inline void set(float &a, const FF b) {
  a = b.x;
}

// Double single assignment
  inline void set(FF &a, const FF b) {
  a = b;
}

// Float assignment
  inline void set(float &a, const float b) {
  a = b;
}

// This function computes b = -a.
 inline FF operator-(const FF a) {
  FF b;
  b.x = -a.x;
  b.y = -a.y;
  return b;
}

// Based on dsadd from DSFUN90, analysis by Norbert Juffa from NVIDIA.
// For `a` and `b` of opposite sign whose magnitude is within a factor of two
// of each other either variant below loses accuracy. Otherwise the result
// is within 1.5ULPs of the correctly rounded result with 48-bit mantissa.
// This function computes c = a + b.
 inline FF operator+(const FF a, const FF b)
{
  FF c;
  float t1, e, t2;
  // Compute dsa + dsb using Knuth's trick.
  t1 = a.x + b.x;
  e = t1 - a.x;
  t2 = ((b.x - e) + (a.x - (t1 - e))) + a.y + b.y;
  // The result is t1 + t2, after normalization.
  c.x = e = t1 + t2;
  c.y = t2 - (e - t1);
  return c;
}

/**
 * Based on dssub from DSFUN90
 * This function computes c = a - b.
 */
 inline
FF operator-(const FF a, const FF b)
{
  FF c;
  float t1, e, t2;
  // Compute dsa - dsb using Knuth's trick.
  t1 = a.x - b.x;
  e = t1 - a.x;
  t2 = ((-b.x - e) + (a.x - (t1 - e))) + a.y - b.y;
  // The result is t1 + t2, after normalization.
  c.x = e = t1 + t2;
  c.y = t2 - (e - t1);
  return c;
}

/**
 * This function multiplies DS numbers A and B to yield the DS product C.
 * Based on: Guillaume Da Graça, David Defour. Implementation of Float-Float
 * Operators on Graphics Hardware. RNC'7, pp. 23-32, 2006.
 *
 */
 inline
FF operator*(const FF a, const FF b)
{
  FF c;
  float up, u1, u2, v1, v2, mh, ml;
  uint32_t tmp;

  // This splits a.x and b.x into high-order and low-order words.
  tmp = (*(uint32_t *)&a.x) & ~0xFFF; // Bit-style splitting from Reimar
  u1 = *(float *)&tmp;
  u2  = a.x - u1;
  tmp = (*(uint32_t *)&b.x) & ~0xFFF;
  v1 = *(float *)&tmp;
  v2  = b.x - v1;
  // Multilply a.x * b.x using Dekker's method.
  mh  = multiply(a.x, b.x);
  ml  = (((multiply(u1, v1) - mh) + multiply(u1, v2))
      + multiply(u2, v1))      + multiply(u2, v2);
  // Compute a.x * b.y + a.y * b.x
  ml  = (multiply(a.x, b.y) + multiply(a.y, b.x)) + ml;
  // The result is mh + ml, after normalization.
  c.x = up = mh + ml;
  c.y = (mh - up) + ml;
  return c;
}

/**
 *  Based on dsdiv from DSFUN90.
 *  This function divides the DS number A by the DS number B to yield the DS
 *  quotient DSC.
 */
 inline FF operator/(const FF a, const FF b)
{
  FF c;
  float s1, cona, conb, a1, b1, a2, b2, c11, c21;
  float c2, t1, e, t2, t12, t22, t11, t21, s2;

  // Compute a DP approximation to the quotient.
  s1 = a.x / b.x;
  // This splits s1 and b.x into high-order and low-order words.
  cona = multiply(s1, 8193.0f);
  conb = multiply(b.x, 8193.0f);
  a1 = cona - (cona - s1);
  b1 = conb - (conb - b.x);
  a2 = s1 - a1;
  b2 = b.x - b1;
  // Multiply s1 * dsb(1) using Dekker's method.
  c11 = multiply(s1, b.x);
  c21 = (((multiply(a1, b1) - c11) + multiply(a1, b2))
      + multiply(a2, b1))       + multiply(a2, b2);
  // Compute s1 * b.y (only high-order word is needed).
  c2 = s1 * b.y;
  // Compute (c11, c21) + c2 using Knuth's trick.
  t1 = c11 + c2;
  e = t1 - c11;
  t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;
  // The result is t1 + t2, after normalization.
  t12 = t1 + t2;
  t22 = t2 - (t12 - t1);
  // Compute dsa - (t12, t22) using Knuth's trick.
  t11 = a.x - t12;
  e = t11 - a.x;
  t21 = ((-t12 - e) + (a.x - (t11 - e))) + a.y - t22;
  // Compute high-order word of (t11, t21) and divide by b.x.
  s2 = (t11 + t21) / b.x;
  // The result is s1 + s2, after normalization.
  c.x = s1 + s2;
  c.y = s2 - (c.x - s1);
  return c;
}



#endif // _FLOAT2_H_
