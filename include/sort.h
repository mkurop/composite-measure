//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sort.h
//
//  Code generation for function 'sort'
//


#ifndef SORT_H
#define SORT_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void sort(coder::array<double, 2U> &x);

#endif

// End of code generation (sort.h)
