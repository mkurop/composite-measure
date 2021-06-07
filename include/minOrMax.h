//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  minOrMax.h
//
//  Code generation for function 'minOrMax'
//


#ifndef MINORMAX_H
#define MINORMAX_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern double b_maximum(const coder::array<double, 2U> &x);
extern void c_maximum(const double x_data[], const int x_size[2], double *ex,
                      int *idx);
extern double maximum(const double x[25]);

#endif

// End of code generation (minOrMax.h)
