#ifndef _F2_BASIC_H_
#define _F2_BASIC_H_

#ifndef _XGPU_
  #define __device__
  #define __host__
  float __fmul_rn(float a, float b) { return a * b;}
#endif



/**
 * Debug report function, for tracing conversions taking place.
 */
//void report(const char* what) {
//  cout << what << endl;
//}

/**
 * Float2's are:
 *   explicitly constructed from `double`,
 *   have `operator=`,
 *   have `operator double`
 */
struct Float2 {
  float x;
  float y;

__host__ __device__    
  Float2() : x(0), y(0) {}

__host__ __device__    
  explicit Float2(double d) {
   // report("creation from double");
    x = (float)d;
    y = (float)(d - x);
  }

__host__ __device__    
  Float2& operator=(const Float2& that)  {
    //report("operator=");
    x = that.x;
    y = that.y;
    return *this;
  }

__host__ __device__
  operator double() const {
    //report("conversion f2->double");
    return (double)x + (double)y;
  }
};

#endif // _F2_BASIC_H_
