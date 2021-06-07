//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  filter.cpp
//
//  Code generation for function 'filter'
//


// Include files
#include <filter.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <gencoswin.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>

// Function Definitions
void filter(const coder::array<double, 2U> &b, const coder::array<double, 2U> &a,
            const coder::array<double, 2U> &x, coder::array<double, 2U> &y)
{
  coder::array<double, 1U> b_b;
  int loop_ub;
  int naxpy;
  coder::array<double, 2U> c_b;
  coder::array<double, 2U> b_a;
  int na;
  int nb;
  double a1;
  int k;
  int nx;
  coder::array<double, 1U> b_y1;
  b_b.set_size(x.size(1));
  loop_ub = x.size(1);
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    b_b[naxpy] = x[naxpy];
  }

  c_b.set_size(1, b.size(1));
  loop_ub = b.size(0) * b.size(1);
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    c_b[naxpy] = b[naxpy];
  }

  b_a.set_size(1, a.size(1));
  loop_ub = a.size(0) * a.size(1);
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    b_a[naxpy] = a[naxpy];
  }

  na = a.size(1);
  nb = b.size(1);
  a1 = a[0];
  if ((!rtIsInf(a[0])) && (!rtIsNaN(a[0])) && (!(a[0] == 0.0)) && (a[0] != 1.0))
  {
    for (k = 0; k < nb; k++) {
      c_b[k] = c_b[k] / a1;
    }

    for (k = 2; k <= na; k++) {
      b_a[k - 1] = b_a[k - 1] / a1;
    }

    b_a[0] = 1.0;
  }

  nx = b_b.size(0) - 1;
  loop_ub = b_b.size(0);
  b_y1.set_size(b_b.size(0));
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    b_y1[naxpy] = 0.0;
  }

  for (k = 0; k <= nx; k++) {
    int j;
    int y1_tmp;
    loop_ub = nx - k;
    naxpy = loop_ub + 1;
    if (naxpy >= nb) {
      naxpy = nb;
    }

    for (j = 0; j < naxpy; j++) {
      y1_tmp = k + j;
      b_y1[y1_tmp] = b_y1[y1_tmp] + b_b[k] * c_b[j];
    }

    naxpy = na - 1;
    if (loop_ub < naxpy) {
      naxpy = loop_ub;
    }

    a1 = -b_y1[k];
    for (j = 0; j < naxpy; j++) {
      y1_tmp = (k + j) + 1;
      b_y1[y1_tmp] = b_y1[y1_tmp] + a1 * b_a[j + 1];
    }
  }

  y.set_size(1, b_y1.size(0));
  loop_ub = b_y1.size(0);
  for (naxpy = 0; naxpy < loop_ub; naxpy++) {
    y[naxpy] = b_y1[naxpy];
  }
}

// End of code generation (filter.cpp)
