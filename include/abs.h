//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  abs.h
//
//  Code generation for function 'abs'
//


#ifndef ABS_H
#define ABS_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void b_abs(const creal_T x_data[], const int x_size[2], double y_data[],
                  int y_size[2]);

#endif

// End of code generation (abs.h)
