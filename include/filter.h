//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  filter.h
//
//  Code generation for function 'filter'
//


#ifndef FILTER_H
#define FILTER_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void filter(const coder::array<double, 2U> &b, const coder::array<double,
                   2U> &a, const coder::array<double, 2U> &x, coder::array<
                   double, 2U> &y);

#endif

// End of code generation (filter.h)
