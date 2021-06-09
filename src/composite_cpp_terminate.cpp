
// Include files
#include <composite_cpp_terminate.h>
#include <composite_cpp.h>
#include <composite_cpp_data.h>
#include <rt_nonfinite.h>

// Function Definitions
void composite_cpp_terminate()
{
  omp_destroy_nest_lock(&emlrtNestLockGlobal);
  isInitialized_composite_cpp = false;
}

// End of code generation (composite_cpp_terminate.cpp)
