//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  log2.h
//
//  Code generation for function 'log2'
//


#ifndef LOG2_H
#define LOG2_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern double b_log2(double x);

#endif

// End of code generation (log2.h)
