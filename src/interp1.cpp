//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  interp1.cpp
//
//  Code generation for function 'interp1'
//


// Include files
#include <interp1.h>
#include <composite_cpp.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <cstring>

// Function Definitions
double interp1(const double varargin_1[26], const double varargin_2[26])
{
  double Vq;
  double y[26];
  double x[26];
  int low_i;
  std::memcpy(&y[0], &varargin_2[0], 26U * sizeof(double));
  std::memcpy(&x[0], &varargin_1[0], 26U * sizeof(double));
  low_i = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (low_i < 26) {
      if (rtIsNaN(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      double xtmp;
      if (varargin_1[1] < varargin_1[0]) {
        for (low_i = 0; low_i < 13; low_i++) {
          xtmp = x[low_i];
          x[low_i] = x[25 - low_i];
          x[25 - low_i] = xtmp;
          xtmp = y[low_i];
          y[low_i] = y[25 - low_i];
          y[25 - low_i] = xtmp;
        }
      }

      Vq = rtNaN;
      if ((!(1000.0 > x[25])) && (!(1000.0 < x[0]))) {
        int low_ip1;
        int high_i;
        low_i = 1;
        low_ip1 = 2;
        high_i = 26;
        while (high_i > low_ip1) {
          int mid_i;
          mid_i = (low_i + high_i) >> 1;
          if (1000.0 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        xtmp = x[low_i - 1];
        xtmp = (1000.0 - xtmp) / (x[low_i] - xtmp);
        if (xtmp == 0.0) {
          Vq = y[low_i - 1];
        } else if (xtmp == 1.0) {
          Vq = y[low_i];
        } else if (y[low_i - 1] == y[low_i]) {
          Vq = y[low_i - 1];
        } else {
          Vq = (1.0 - xtmp) * y[low_i - 1] + xtmp * y[low_i];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return Vq;
}

void interp1Linear(const double y[26], const coder::array<double, 2U> &xi, coder::
                   array<double, 2U> &yi, const double varargin_1[26])
{
  double minx;
  double maxx;
  int ub_loop;
  double r_tmp;
  int low_i;
  int low_ip1;
  int high_i;
  int mid_i;
  double r;
  minx = varargin_1[0];
  maxx = varargin_1[25];
  ub_loop = xi.size(1) - 1;

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(r_tmp,low_i,low_ip1,high_i,mid_i,r)

  for (int k = 0; k <= ub_loop; k++) {
    r_tmp = xi[k];
    if (rtIsNaN(r_tmp)) {
      yi[k] = rtNaN;
    } else {
      if ((!(r_tmp > maxx)) && (!(r_tmp < minx))) {
        low_i = 1;
        low_ip1 = 2;
        high_i = 26;
        while (high_i > low_ip1) {
          mid_i = (low_i + high_i) >> 1;
          if (xi[k] >= varargin_1[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }

        r_tmp = varargin_1[low_i - 1];
        r = (xi[k] - r_tmp) / (varargin_1[low_i] - r_tmp);
        if (r == 0.0) {
          yi[k] = y[low_i - 1];
        } else if (r == 1.0) {
          yi[k] = y[low_i];
        } else {
          r_tmp = y[low_i - 1];
          if (r_tmp == y[low_i]) {
            yi[k] = r_tmp;
          } else {
            yi[k] = (1.0 - r) * r_tmp + r * y[low_i];
          }
        }
      }
    }
  }
}

// End of code generation (interp1.cpp)
