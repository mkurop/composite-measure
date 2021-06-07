//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  composite_cpp.h
//
//  Code generation for function 'composite_cpp'
//


#ifndef COMPOSITE_CPP_H
#define COMPOSITE_CPP_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void composite_cpp(const coder::array<double, 1U> &cleanFile, const coder::
  array<double, 1U> &enhancedFile, double Srate1, double Srate2, double
  Csig_data[], int Csig_size[2], double Cbak_data[], int Cbak_size[2], double
  Covl_data[], int Covl_size[2]);

#endif

// End of code generation (composite_cpp.h)
