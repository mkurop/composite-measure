//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  find.cpp
//
//  Code generation for function 'find'
//


// Include files
#include <find.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>

// Function Definitions
void eml_find(const coder::array<boolean_T, 2U> &x, coder::array<int, 2U> &i)
{
  int nx;
  int idx;
  int ii;
  boolean_T exitg1;
  nx = x.size(1);
  idx = 0;
  i.set_size(1, x.size(1));
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x[ii]) {
      idx++;
      i[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }

  if (1 > idx) {
    idx = 0;
  }

  i.set_size(i.size(0), idx);
}

// End of code generation (find.cpp)
