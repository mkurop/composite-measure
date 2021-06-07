//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  hann.h
//
//  Code generation for function 'hann'
//


#ifndef HANN_H
#define HANN_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void hann(double varargin_1, double w_data[], int w_size[1]);

#endif

// End of code generation (hann.h)
