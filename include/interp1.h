//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  interp1.h
//
//  Code generation for function 'interp1'
//


#ifndef INTERP1_H
#define INTERP1_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern double interp1(const double varargin_1[26], const double varargin_2[26]);
extern void interp1Linear(const double y[26], const coder::array<double, 2U> &xi,
  coder::array<double, 2U> &yi, const double varargin_1[26]);

#endif

// End of code generation (interp1.h)
