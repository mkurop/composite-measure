//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sortIdx.h
//
//  Code generation for function 'sortIdx'
//


#ifndef SORTIDX_H
#define SORTIDX_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void merge_block(coder::array<int, 2U> &idx, coder::array<double, 2U> &x,
  int offset, int n, int preSortLevel, coder::array<int, 1U> &iwork, coder::
  array<double, 1U> &xwork);

#endif

// End of code generation (sortIdx.h)
