//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  pesq_original_cpp.h
//
//  Code generation for function 'pesq_original_cpp'
//


#ifndef PESQ_ORIGINAL_CPP_H
#define PESQ_ORIGINAL_CPP_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void pesq_original_cpp(const coder::array<double, 1U> &ref_data, double
  ref_sampling_rate, const coder::array<double, 1U> &deg_data, double
  scores_data[], int scores_size[2]);

#endif

// End of code generation (pesq_original_cpp.h)
