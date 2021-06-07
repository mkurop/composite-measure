//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  gencoswin.cpp
//
//  Code generation for function 'gencoswin'
//


// Include files
#include <gencoswin.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <cmath>

// Function Definitions
void calc_window(double m, double n, double w_data[], int w_size[1])
{
  coder::array<double, 2U> y;
  int nx;
  coder::array<double, 1U> b_y;
  int k;
  if (rtIsNaN(m - 1.0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else if (m - 1.0 < 0.0) {
    y.set_size(1, 0);
  } else if (rtIsInf(m - 1.0) && (0.0 == m - 1.0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else {
    nx = static_cast<int>(std::floor(m - 1.0));
    y.set_size(1, (nx + 1));
    for (k = 0; k <= nx; k++) {
      y[k] = k;
    }
  }

  b_y.set_size(y.size(1));
  nx = y.size(1);
  for (k = 0; k < nx; k++) {
    b_y[k] = y[k];
  }

  nx = b_y.size(0);
  for (k = 0; k < nx; k++) {
    b_y[k] = 6.2831853071795862 * (b_y[k] / (n - 1.0));
  }

  nx = b_y.size(0);
  for (k = 0; k < nx; k++) {
    b_y[k] = std::cos(b_y[k]);
  }

  w_size[0] = b_y.size(0);
  nx = b_y.size(0);
  for (k = 0; k < nx; k++) {
    w_data[k] = 0.5 - 0.5 * b_y[k];
  }
}

// End of code generation (gencoswin.cpp)
