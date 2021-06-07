//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  ifft.cpp
//
//  Code generation for function 'ifft'
//


// Include files
#include <ifft.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <fft.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <cmath>

// Function Definitions
void ifft(const coder::array<creal_T, 2U> &x, double varargin_1, coder::array<
          creal_T, 2U> &y)
{
  int nfft_tmp;
  int N2blue;
  int nd2;
  coder::array<double, 2U> costab1q;
  coder::array<double, 2U> costab;
  coder::array<double, 2U> sintab;
  coder::array<double, 2U> sintabinv;
  coder::array<creal_T, 1U> wwc;
  coder::array<creal_T, 1U> yCol;
  int rt;
  coder::array<creal_T, 1U> fv;
  coder::array<creal_T, 1U> b_fv;
  nfft_tmp = static_cast<int>(varargin_1);
  if ((x.size(1) == 0) || (0 == nfft_tmp)) {
    y.set_size(1, nfft_tmp);
    for (int nInt2m1 = 0; nInt2m1 < nfft_tmp; nInt2m1++) {
      y[nInt2m1].re = 0.0;
      y[nInt2m1].im = 0.0;
    }
  } else {
    int useRadix2_tmp;
    boolean_T useRadix2;
    int nInt2m1;
    double nt_im;
    int xidx;
    int k;
    useRadix2_tmp = nfft_tmp - 1;
    useRadix2 = ((nfft_tmp > 0) && ((nfft_tmp & useRadix2_tmp) == 0));
    FFTImplementationCallback::get_algo_sizes((static_cast<int>(varargin_1)),
      (useRadix2), (&N2blue), (&nd2));
    nt_im = 6.2831853071795862 / static_cast<double>(nd2);
    xidx = nd2 / 2 / 2;
    costab1q.set_size(1, (xidx + 1));
    costab1q[0] = 1.0;
    nd2 = xidx / 2 - 1;
    for (k = 0; k <= nd2; k++) {
      costab1q[k + 1] = std::cos(nt_im * (static_cast<double>(k) + 1.0));
    }

    nInt2m1 = nd2 + 2;
    nd2 = xidx - 1;
    for (k = nInt2m1; k <= nd2; k++) {
      costab1q[k] = std::sin(nt_im * static_cast<double>(xidx - k));
    }

    costab1q[xidx] = 0.0;
    if (!useRadix2) {
      xidx = costab1q.size(1) - 1;
      nd2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, (nd2 + 1));
      sintab.set_size(1, (nd2 + 1));
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, (nd2 + 1));
      for (k = 0; k < xidx; k++) {
        sintabinv[k + 1] = costab1q[(xidx - k) - 1];
      }

      nInt2m1 = costab1q.size(1);
      for (k = nInt2m1; k <= nd2; k++) {
        sintabinv[k] = costab1q[k - xidx];
      }

      for (k = 0; k < xidx; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(xidx - k) - 1];
      }

      nInt2m1 = costab1q.size(1);
      for (k = nInt2m1; k <= nd2; k++) {
        costab[k] = -costab1q[nd2 - k];
        sintab[k] = -costab1q[k - xidx];
      }
    } else {
      xidx = costab1q.size(1) - 1;
      nd2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, (nd2 + 1));
      sintab.set_size(1, (nd2 + 1));
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (k = 0; k < xidx; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = costab1q[(xidx - k) - 1];
      }

      nInt2m1 = costab1q.size(1);
      for (k = nInt2m1; k <= nd2; k++) {
        costab[k] = -costab1q[nd2 - k];
        sintab[k] = costab1q[k - xidx];
      }

      sintabinv.set_size(1, 0);
    }

    if (useRadix2) {
      nd2 = x.size(1);
      wwc = x.reshape(nd2);
      FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (nfft_tmp), (costab),
        (sintab), (yCol));
      if (yCol.size(0) > 1) {
        nt_im = 1.0 / static_cast<double>(yCol.size(0));
        nd2 = yCol.size(0);
        for (nInt2m1 = 0; nInt2m1 < nd2; nInt2m1++) {
          yCol[nInt2m1].re = nt_im * yCol[nInt2m1].re;
          yCol[nInt2m1].im = nt_im * yCol[nInt2m1].im;
        }
      }
    } else {
      int idx;
      double nt_re;
      nInt2m1 = (nfft_tmp + nfft_tmp) - 1;
      wwc.set_size(nInt2m1);
      idx = nfft_tmp;
      rt = 0;
      wwc[useRadix2_tmp].re = 1.0;
      wwc[useRadix2_tmp].im = 0.0;
      nd2 = nfft_tmp << 1;
      for (k = 0; k <= nfft_tmp - 2; k++) {
        xidx = ((k + 1) << 1) - 1;
        if (nd2 - rt <= xidx) {
          rt += xidx - nd2;
        } else {
          rt += xidx;
        }

        nt_im = 3.1415926535897931 * static_cast<double>(rt) / static_cast<
          double>(nfft_tmp);
        if (nt_im == 0.0) {
          nt_re = 1.0;
          nt_im = 0.0;
        } else {
          nt_re = std::cos(nt_im);
          nt_im = std::sin(nt_im);
        }

        wwc[idx - 2].re = nt_re;
        wwc[idx - 2].im = -nt_im;
        idx--;
      }

      idx = 0;
      nInt2m1--;
      for (k = nInt2m1; k >= nfft_tmp; k--) {
        wwc[k] = wwc[idx];
        idx++;
      }

      yCol.set_size(nfft_tmp);
      if (nfft_tmp > x.size(1)) {
        yCol.set_size(nfft_tmp);
        for (nInt2m1 = 0; nInt2m1 < nfft_tmp; nInt2m1++) {
          yCol[nInt2m1].re = 0.0;
          yCol[nInt2m1].im = 0.0;
        }
      }

      nd2 = x.size(1);
      if (nfft_tmp < nd2) {
        nd2 = nfft_tmp;
      }

      xidx = 0;
      for (k = 0; k < nd2; k++) {
        nInt2m1 = (nfft_tmp + k) - 1;
        yCol[k].re = wwc[nInt2m1].re * x[xidx].re + wwc[nInt2m1].im * x[xidx].im;
        yCol[k].im = wwc[nInt2m1].re * x[xidx].im - wwc[nInt2m1].im * x[xidx].re;
        xidx++;
      }

      nInt2m1 = nd2 + 1;
      for (k = nInt2m1; k <= nfft_tmp; k++) {
        yCol[k - 1].re = 0.0;
        yCol[k - 1].im = 0.0;
      }

      FFTImplementationCallback::r2br_r2dit_trig_impl((yCol), (N2blue), (costab),
        (sintab), (fv));
      FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (N2blue), (costab),
        (sintab), (b_fv));
      b_fv.set_size(fv.size(0));
      nd2 = fv.size(0);
      for (nInt2m1 = 0; nInt2m1 < nd2; nInt2m1++) {
        nt_im = fv[nInt2m1].re * b_fv[nInt2m1].im + fv[nInt2m1].im *
          b_fv[nInt2m1].re;
        b_fv[nInt2m1].re = fv[nInt2m1].re * b_fv[nInt2m1].re - fv[nInt2m1].im *
          b_fv[nInt2m1].im;
        b_fv[nInt2m1].im = nt_im;
      }

      FFTImplementationCallback::r2br_r2dit_trig_impl((b_fv), (N2blue), (costab),
        (sintabinv), (fv));
      if (fv.size(0) > 1) {
        nt_im = 1.0 / static_cast<double>(fv.size(0));
        nd2 = fv.size(0);
        for (nInt2m1 = 0; nInt2m1 < nd2; nInt2m1++) {
          fv[nInt2m1].re = nt_im * fv[nInt2m1].re;
          fv[nInt2m1].im = nt_im * fv[nInt2m1].im;
        }
      }

      idx = 0;
      nt_re = static_cast<int>(varargin_1);
      nInt2m1 = wwc.size(0);
      for (k = nfft_tmp; k <= nInt2m1; k++) {
        double re;
        yCol[idx].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].
          im;
        yCol[idx].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].
          re;
        if (yCol[idx].im == 0.0) {
          re = yCol[idx].re / nt_re;
          nt_im = 0.0;
        } else if (yCol[idx].re == 0.0) {
          re = 0.0;
          nt_im = yCol[idx].im / nt_re;
        } else {
          re = yCol[idx].re / nt_re;
          nt_im = yCol[idx].im / nt_re;
        }

        yCol[idx].re = re;
        yCol[idx].im = nt_im;
        idx++;
      }
    }

    y.set_size(1, nfft_tmp);
    for (nInt2m1 = 0; nInt2m1 < nfft_tmp; nInt2m1++) {
      y[nInt2m1] = yCol[nInt2m1];
    }
  }
}

// End of code generation (ifft.cpp)
