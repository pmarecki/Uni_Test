/**
 *  @brief Stopers for various architectures. 
 *  
 *  Defined: 
 *   \li PosixStoper        -- (Linux/MacOS) 
 *   \li CudaStoper(CUDA)   -- synchronizes with GPU
 *   \li NanoStoper (Linux) -- best resolution; 
 *                          clock_gettime(CLOCK_THREAD_CPUTIME_ID,..)
 *   \li NanoStoper (MacOS) -- mach_absolute_time()
 *
 *  @note For CUDA user must #define XGPU
 *  @note On Linux compile with `-lrt`
 *
 *  Version 1.0
 *
 */

 
//    Usage (for all stopers):
//      PosixStoper xx;
//      /* some work */
//      xx.Snap();
//      //cout << xx.LastDt();
//      xx.Print();

#define MAX_NUMBER_SNAPS 10000 //maximum no. of time-measurements per stoper


#ifndef STOPPER_H_
#define STOPPER_H_

#include <sys/time.h>   //time functions, gettimeofday()
#include <inttypes.h>   //int64_t  
#include <cstdio>
#include <cassert> 
#include <string.h>     //memset
#include <time.h>       //clock_gettime

		
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif



/////////////////////////////
// General stoper class.
//
// Members: Snap(), Diff(int t2, int t1), must be implemented in true stopers.

class StoperBase {        
public:    
    
  StoperBase() : position_(0) { memcpy(unit_name_,&"usec",5);}
    virtual ~StoperBase() {}

    const char* unit_name() {return &unit_name_[0];}
    void Reset()  { position_ = 0;}
    virtual void Snap() { position_++;} //records current time  
    virtual float Diff(int t2, int t1) { return t2 - t1;}
  
    virtual float LastDt() {  //Must have measured >=2 values.
    return (position_ < 2) ?  0 : Diff(position_ - 1, position_ - 2); 
  }
    void Print(void) {
    if (position_ == 0) return;
    for(int i=1; i < position_; i++) 
      printf("[%i:%i]\t=%6.3f %s\n", i-1, i, Diff(i, i - 1), unit_name_);
  }
protected:
  int position_;
  char unit_name_[5];
};



/////////////////////////////
// Implementation part


/////////////////////////////
// POSIX Stoper 
// MacOS Latency = 100ns, +/- 0.5ns; granularity = 1us
// Linux Latency = 1.3us, +/- 40ns 
// Uses `gettimeofday()`.
class PosixStoper : public StoperBase {
 public:
  PosixStoper() { Snap();}
  virtual inline void Snap()  {
    gettimeofday(&tv_[position_++], NULL);
  }
  virtual float Diff(int t2, int t1) {
    float elapsed = (tv_[t2].tv_sec - tv_[t1].tv_sec) * 1000000.0 +
        (tv_[t2].tv_usec-tv_[t1].tv_usec);
    return elapsed;
  }

 private:
  timeval tv_[MAX_NUMBER_SNAPS];
};



/////////////////////////////
// Nanosecond POSIX stoper.
// Uses `clock_gettime()`. Latency = 400ns, Accuracy = 40ns
#ifndef __MACH__ //not defined on Mach
class NanoStoper : public StoperBase {
public:
	NanoStoper()  { Snap();}
  inline virtual void Snap() { 
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_[position_++]);
  }
  inline virtual float Diff(int t2, int t1) {
		float elapsed = ((tv_[t2].tv_sec - tv_[t1].tv_sec) * 1000000.0 + 
		                 (tv_[t2].tv_nsec- tv_[t1].tv_nsec) / 1000.0);
		return elapsed;		
	}
 private:
  timespec tv_[MAX_NUMBER_SNAPS];
};
#endif // not __MACH__



/////////////////////////////
// Nanosecod MacOS (Mach) stoper.
// Uses `mach_absolute_time()`. Latency = 65ns, Accuracy = 2ns. 
#ifdef __MACH__
class NanoStoper : public StoperBase {
 public:
	NanoStoper() {
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    conversion_factor_ = (double)timebase.numer / (double)timebase.denom;
	  Snap();
	}
  ~NanoStoper() {}
  virtual inline void Snap()   { tv_[position_++] = mach_absolute_time();}
  
  inline virtual float Diff(int t2, int t1) {
		double elapsed = (double)(tv_[t2] - tv_[t1]) * conversion_factor_ /1000.0;
		return elapsed;		
	}
 private:
  double conversion_factor_; 
  uint64_t tv_[MAX_NUMBER_SNAPS];
};
#endif //__MACH__


/////////////////////////////
// GPU-CUDA stoper to be placed in host (CPU) code.
// Latency = 16us, +/- 1.5us
// Invoked by host, creates `cudeEvent_t` on device, and synchronizes with them.
#ifdef XGPU
class CudaStoper : public StoperBase {  
 public:
  CudaStoper()  { Snap();}
  // synchronizes device as if cudaDeviceSynchronize().
   inline virtual void Snap() {   
    cudaEventCreate(&tv_[position_]);
    cudaEventRecord(tv_[position_], 0);
    cudaEventSynchronize(tv_[position_++]);
  }
   inline virtual float Diff(int t2, int t1) {
    float elapsed;
    cudaEventElapsedTime(&elapsed, tv_[t1], tv_[t2]); //elapsed in ms
		return elapsed * 1000.0;		
	}
 private:
  cudaEvent_t tv_[MAX_NUMBER_SNAPS];
};
#endif //XGPU




#endif //STOPPER_H_


