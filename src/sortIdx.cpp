//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sortIdx.cpp
//
//  Code generation for function 'sortIdx'
//


// Include files
#include <sortIdx.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <gencoswin.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>

// Function Declarations
static void merge(coder::array<int, 2U> &idx, coder::array<double, 2U> &x, int
                  offset, int np, int nq, coder::array<int, 1U> &iwork, coder::
                  array<double, 1U> &xwork);

// Function Definitions
static void merge(coder::array<int, 2U> &idx, coder::array<double, 2U> &x, int
                  offset, int np, int nq, coder::array<int, 1U> &iwork, coder::
                  array<double, 1U> &xwork)
{
  if (nq != 0) {
    int n_tmp;
    int j;
    int p;
    int iout;
    int q;
    n_tmp = np + nq;
    for (j = 0; j < n_tmp; j++) {
      iout = offset + j;
      iwork[j] = idx[iout];
      xwork[j] = x[iout];
    }

    p = 0;
    q = np;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[q]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[q];
        x[iout] = xwork[q];
        if (q + 1 < n_tmp) {
          q++;
        } else {
          q = iout - p;
          for (j = p + 1; j <= np; j++) {
            iout = q + j;
            idx[iout] = iwork[j - 1];
            x[iout] = xwork[j - 1];
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

void merge_block(coder::array<int, 2U> &idx, coder::array<double, 2U> &x, int
                 offset, int n, int preSortLevel, coder::array<int, 1U> &iwork,
                 coder::array<double, 1U> &xwork)
{
  int nPairs;
  int bLen;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    int tailOffset;
    int nTail;
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

// End of code generation (sortIdx.cpp)
