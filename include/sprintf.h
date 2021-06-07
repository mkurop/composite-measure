//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sprintf.h
//
//  Code generation for function 'sprintf'
//


#ifndef SPRINTF_H
#define SPRINTF_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void b_sprintf(const char varargin_1_data[], const int varargin_1_size[2],
                      coder::array<char, 2U> &str);

#endif

// End of code generation (sprintf.h)
