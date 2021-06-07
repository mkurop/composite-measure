//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  relop.h
//
//  Code generation for function 'relop'
//


#ifndef RELOP_H
#define RELOP_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "composite_cpp_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern boolean_T iseq(double x, double y);

#endif

// End of code generation (relop.h)
