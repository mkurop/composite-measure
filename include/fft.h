//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  fft.h
//
//  Code generation for function 'fft'
//


#ifndef FFT_H
#define FFT_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void b_fft(const coder::array<double, 2U> &x, double varargin_1, coder::
                  array<creal_T, 2U> &y);
extern void c_fft(const double x_data[], const int x_size[2], creal_T y_data[],
                  int y_size[2]);
extern void fft(const coder::array<double, 1U> &x, double varargin_1, coder::
                array<creal_T, 1U> &y);

#endif

// End of code generation (fft.h)
