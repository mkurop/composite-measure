//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  abs.cpp
//
//  Code generation for function 'abs'
//


// Include files
#include <abs.h>
#include <composite_cpp.h>
#include <composite_cpp_rtwutil.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>

// Function Definitions
void b_abs(const creal_T x_data[], const int x_size[2], double y_data[], int
           y_size[2])
{
  int nx;
  nx = x_size[1];
  y_size[0] = 1;
  y_size[1] = static_cast<short>(x_size[1]);
  for (int k = 0; k < nx; k++) {
    y_data[k] = rt_hypotd_snf(x_data[k].re, x_data[k].im);
  }
}

// End of code generation (abs.cpp)
