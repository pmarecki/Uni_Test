/**
 *  Error handling/reporting for CUDA.
 *   Methods/macros:
 *    -- _try(-cuda_command-)
 *    -- CUDA_CHK
 *
 *  Version: 0.8.7b
 */

// Common errors:

//Too Many Resources Requested for Launch -
//  This error means that the number of registers available on multiprocessor
//  is being exceeded. Reduce threads per block to resolve the problem.
//Unspecified launch failure -
//  Translate it to "segmentation fault" for the host code.

#ifndef CUDA_ERROR_
#define CUDA_ERROR_

#include <stdio.h>


/**
 * Output CUDA error-strings.
 *
 * Usage: _try(cudaMalloc(...)); 
 *
 */
#define _try(err)  __checkCudaErrors (err, __FILE__, __LINE__)
inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
  if (cudaSuccess != err) {
    fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n", file, line, 
            (int)err, cudaGetErrorString(err) );
    exit(-1);        
  }
}


/**
 *  Printing `cudaGetLastError()`.
 *  Usage: 
 *        kernel<<<1,1>>(args);
 *        CUDA_CHK
 */
#define CUDA_CHK  __getLastCudaError(__FILE__, __LINE__)
inline void __getLastCudaError(const char *file_name, const int line)
{
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "(file: %s line: %i) cuERR:  (%d) %s.\n", file_name, line,
            (int)err, cudaGetErrorString(err));
    exit(-1);
  }
}


#endif //CUDA_ERROR_

