//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  composite_cpp_initialize.cpp
//
//  Code generation for function 'composite_cpp_initialize'
//


// Include files
#include <composite_cpp_initialize.h>
#include <composite_cpp.h>
#include <composite_cpp_data.h>
#include <rt_nonfinite.h>

// Function Definitions
void composite_cpp_initialize()
{
  rt_InitInfAndNaN();
  omp_init_nest_lock(&emlrtNestLockGlobal);
  isInitialized_composite_cpp = true;
}

// End of code generation (composite_cpp_initialize.cpp)
