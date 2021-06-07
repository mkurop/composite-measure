//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  hann.cpp
//
//  Code generation for function 'hann'
//


// Include files
#include <hann.h>
#include <composite_cpp.h>
#include <composite_cpp_rtwutil.h>
#include <gencoswin.h>
#include <rt_nonfinite.h>
#include <cmath>
#include <cstring>

// Function Definitions
void hann(double varargin_1, double w_data[], int w_size[1])
{
  double L;
  double r;
  double tmp_data[257];
  int tmp_size[1];
  double b_w_data[514];
  if (varargin_1 == std::floor(varargin_1)) {
    L = varargin_1;
  } else {
    L = rt_roundd_snf(varargin_1);
  }

  if (rtIsNaN(L + 1.0) || rtIsInf(L + 1.0)) {
    r = rtNaN;
  } else {
    r = std::fmod(L + 1.0, 2.0);
    if (r == 0.0) {
      r = 0.0;
    }
  }

  if (r == 0.0) {
    int loop_ub;
    int i;
    int i1;
    int b_loop_ub;
    int w_size_idx_0;
    calc_window((L + 1.0) / 2.0, L + 1.0, tmp_data, tmp_size);
    loop_ub = tmp_size[0];
    if (0 <= loop_ub - 1) {
      std::memcpy(&w_data[0], &tmp_data[0], loop_ub * sizeof(double));
    }

    if (2 > tmp_size[0]) {
      i = 0;
      i1 = 1;
      b_loop_ub = -1;
    } else {
      i = tmp_size[0] - 1;
      i1 = -1;
      b_loop_ub = 1;
    }

    loop_ub = div_s32_floor(b_loop_ub - i, i1);
    w_size_idx_0 = (tmp_size[0] + loop_ub) + 1;
    b_loop_ub = tmp_size[0];
    if (0 <= b_loop_ub - 1) {
      std::memcpy(&b_w_data[0], &w_data[0], b_loop_ub * sizeof(double));
    }

    for (b_loop_ub = 0; b_loop_ub <= loop_ub; b_loop_ub++) {
      b_w_data[b_loop_ub + tmp_size[0]] = w_data[i + i1 * b_loop_ub];
    }

    w_size[0] = w_size_idx_0;
    if (0 <= w_size_idx_0 - 1) {
      std::memcpy(&w_data[0], &b_w_data[0], w_size_idx_0 * sizeof(double));
    }
  } else {
    int loop_ub;
    int i;
    int i1;
    int b_loop_ub;
    int w_size_idx_0;
    calc_window(((L + 1.0) + 1.0) / 2.0, L + 1.0, tmp_data, tmp_size);
    loop_ub = tmp_size[0];
    if (0 <= loop_ub - 1) {
      std::memcpy(&w_data[0], &tmp_data[0], loop_ub * sizeof(double));
    }

    if (2 > tmp_size[0] - 1) {
      i = 0;
      i1 = 1;
      b_loop_ub = -1;
    } else {
      i = tmp_size[0] - 2;
      i1 = -1;
      b_loop_ub = 1;
    }

    loop_ub = div_s32_floor(b_loop_ub - i, i1);
    w_size_idx_0 = (tmp_size[0] + loop_ub) + 1;
    b_loop_ub = tmp_size[0];
    if (0 <= b_loop_ub - 1) {
      std::memcpy(&b_w_data[0], &w_data[0], b_loop_ub * sizeof(double));
    }

    for (b_loop_ub = 0; b_loop_ub <= loop_ub; b_loop_ub++) {
      b_w_data[b_loop_ub + tmp_size[0]] = w_data[i + i1 * b_loop_ub];
    }

    w_size[0] = w_size_idx_0;
    if (0 <= w_size_idx_0 - 1) {
      std::memcpy(&w_data[0], &b_w_data[0], w_size_idx_0 * sizeof(double));
    }
  }
}

// End of code generation (hann.cpp)
