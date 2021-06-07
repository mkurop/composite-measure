//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  composite_cpp_rtwutil.h
//
//  Code generation for function 'composite_cpp_rtwutil'
//


#ifndef COMPOSITE_CPP_RTWUTIL_H
#define COMPOSITE_CPP_RTWUTIL_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern int div_s32_floor(int numerator, int denominator);
extern double rt_hypotd_snf(double u0, double u1);
extern double rt_powd_snf(double u0, double u1);
extern double rt_roundd_snf(double u);

#endif

// End of code generation (composite_cpp_rtwutil.h)
