//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sort.cpp
//
//  Code generation for function 'sort'
//


// Include files
#include <sort.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <gencoswin.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <sortIdx.h>

// Function Definitions
void sort(coder::array<double, 2U> &x)
{
  coder::array<int, 2U> idx;
  int i1;
  int ib;
  double x4[4];
  int idx4[4];
  coder::array<int, 1U> iwork;
  coder::array<double, 1U> xwork;
  signed char perm[4];
  int b_iwork[256];
  double b_xwork[256];
  idx.set_size(1, x.size(1));
  i1 = x.size(1);
  for (ib = 0; ib < i1; ib++) {
    idx[ib] = 0;
  }

  if (x.size(1) != 0) {
    int n;
    int b_n;
    int bLen;
    int k;
    int i4;
    int idx_tmp;
    int quartetOffset;
    int i3;
    int nNonNaN;
    n = x.size(1);
    b_n = x.size(1);
    x4[0] = 0.0;
    idx4[0] = 0;
    x4[1] = 0.0;
    idx4[1] = 0;
    x4[2] = 0.0;
    idx4[2] = 0;
    x4[3] = 0.0;
    idx4[3] = 0;
    i1 = x.size(1);
    iwork.set_size(i1);
    for (ib = 0; ib < i1; ib++) {
      iwork[ib] = 0;
    }

    i1 = x.size(1);
    xwork.set_size(i1);
    for (ib = 0; ib < i1; ib++) {
      xwork[ib] = 0.0;
    }

    bLen = 0;
    ib = -1;
    for (k = 0; k < b_n; k++) {
      if (rtIsNaN(x[k])) {
        idx_tmp = (b_n - bLen) - 1;
        idx[idx_tmp] = k + 1;
        xwork[idx_tmp] = x[k];
        bLen++;
      } else {
        ib++;
        idx4[ib] = k + 1;
        x4[ib] = x[k];
        if (ib + 1 == 4) {
          double d;
          double d1;
          quartetOffset = k - bLen;
          if (x4[0] <= x4[1]) {
            i1 = 1;
            ib = 2;
          } else {
            i1 = 2;
            ib = 1;
          }

          if (x4[2] <= x4[3]) {
            i3 = 3;
            i4 = 4;
          } else {
            i3 = 4;
            i4 = 3;
          }

          d = x4[i1 - 1];
          d1 = x4[i3 - 1];
          if (d <= d1) {
            d = x4[ib - 1];
            if (d <= d1) {
              perm[0] = static_cast<signed char>(i1);
              perm[1] = static_cast<signed char>(ib);
              perm[2] = static_cast<signed char>(i3);
              perm[3] = static_cast<signed char>(i4);
            } else if (d <= x4[i4 - 1]) {
              perm[0] = static_cast<signed char>(i1);
              perm[1] = static_cast<signed char>(i3);
              perm[2] = static_cast<signed char>(ib);
              perm[3] = static_cast<signed char>(i4);
            } else {
              perm[0] = static_cast<signed char>(i1);
              perm[1] = static_cast<signed char>(i3);
              perm[2] = static_cast<signed char>(i4);
              perm[3] = static_cast<signed char>(ib);
            }
          } else {
            d1 = x4[i4 - 1];
            if (d <= d1) {
              if (x4[ib - 1] <= d1) {
                perm[0] = static_cast<signed char>(i3);
                perm[1] = static_cast<signed char>(i1);
                perm[2] = static_cast<signed char>(ib);
                perm[3] = static_cast<signed char>(i4);
              } else {
                perm[0] = static_cast<signed char>(i3);
                perm[1] = static_cast<signed char>(i1);
                perm[2] = static_cast<signed char>(i4);
                perm[3] = static_cast<signed char>(ib);
              }
            } else {
              perm[0] = static_cast<signed char>(i3);
              perm[1] = static_cast<signed char>(i4);
              perm[2] = static_cast<signed char>(i1);
              perm[3] = static_cast<signed char>(ib);
            }
          }

          idx_tmp = perm[0] - 1;
          idx[quartetOffset - 3] = idx4[idx_tmp];
          i3 = perm[1] - 1;
          idx[quartetOffset - 2] = idx4[i3];
          ib = perm[2] - 1;
          idx[quartetOffset - 1] = idx4[ib];
          i1 = perm[3] - 1;
          idx[quartetOffset] = idx4[i1];
          x[quartetOffset - 3] = x4[idx_tmp];
          x[quartetOffset - 2] = x4[i3];
          x[quartetOffset - 1] = x4[ib];
          x[quartetOffset] = x4[i1];
          ib = -1;
        }
      }
    }

    i4 = (b_n - bLen) - 1;
    if (ib + 1 > 0) {
      perm[1] = 0;
      perm[2] = 0;
      perm[3] = 0;
      if (ib + 1 == 1) {
        perm[0] = 1;
      } else if (ib + 1 == 2) {
        if (x4[0] <= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] <= x4[1]) {
        if (x4[1] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 3;
        perm[2] = 1;
      } else {
        perm[0] = 3;
        perm[1] = 2;
        perm[2] = 1;
      }

      for (k = 0; k <= ib; k++) {
        idx_tmp = perm[k] - 1;
        i3 = (i4 - ib) + k;
        idx[i3] = idx4[idx_tmp];
        x[i3] = x4[idx_tmp];
      }
    }

    ib = (bLen >> 1) + 1;
    for (k = 0; k <= ib - 2; k++) {
      i1 = (i4 + k) + 1;
      i3 = idx[i1];
      idx_tmp = (b_n - k) - 1;
      idx[i1] = idx[idx_tmp];
      idx[idx_tmp] = i3;
      x[i1] = xwork[idx_tmp];
      x[idx_tmp] = xwork[i1];
    }

    if ((bLen & 1) != 0) {
      ib += i4;
      x[ib] = xwork[ib];
    }

    nNonNaN = n - bLen;
    i1 = 2;
    if (nNonNaN > 1) {
      if (n >= 256) {
        int nBlocks;
        nBlocks = nNonNaN >> 8;
        if (nBlocks > 0) {
          for (int b = 0; b < nBlocks; b++) {
            quartetOffset = (b << 8) - 1;
            for (int b_b = 0; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              b_n = bLen << 1;
              n = 256 >> (b_b + 3);
              for (k = 0; k < n; k++) {
                i3 = (quartetOffset + k * b_n) + 1;
                for (i1 = 0; i1 < b_n; i1++) {
                  ib = i3 + i1;
                  b_iwork[i1] = idx[ib];
                  b_xwork[i1] = x[ib];
                }

                i4 = 0;
                i1 = bLen;
                ib = i3 - 1;
                int exitg1;
                do {
                  exitg1 = 0;
                  ib++;
                  if (b_xwork[i4] <= b_xwork[i1]) {
                    idx[ib] = b_iwork[i4];
                    x[ib] = b_xwork[i4];
                    if (i4 + 1 < bLen) {
                      i4++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    idx[ib] = b_iwork[i1];
                    x[ib] = b_xwork[i1];
                    if (i1 + 1 < b_n) {
                      i1++;
                    } else {
                      ib -= i4;
                      for (i1 = i4 + 1; i1 <= bLen; i1++) {
                        idx_tmp = ib + i1;
                        idx[idx_tmp] = b_iwork[i1 - 1];
                        x[idx_tmp] = b_xwork[i1 - 1];
                      }

                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }

          i1 = nBlocks << 8;
          ib = nNonNaN - i1;
          if (ib > 0) {
            merge_block(idx, x, i1, ib, 2, iwork, xwork);
          }

          i1 = 8;
        }
      }

      merge_block(idx, x, 0, nNonNaN, i1, iwork, xwork);
    }
  }
}

// End of code generation (sort.cpp)
