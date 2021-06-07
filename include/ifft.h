//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  ifft.h
//
//  Code generation for function 'ifft'
//


#ifndef IFFT_H
#define IFFT_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void ifft(const coder::array<creal_T, 2U> &x, double varargin_1, coder::
                 array<creal_T, 2U> &y);

#endif

// End of code generation (ifft.h)
