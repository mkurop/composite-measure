//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  find.h
//
//  Code generation for function 'find'
//


#ifndef FIND_H
#define FIND_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void eml_find(const coder::array<boolean_T, 2U> &x, coder::array<int, 2U>
                     &i);

#endif

// End of code generation (find.h)
