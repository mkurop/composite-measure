//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  FFTImplementationCallback.cpp
//
//  Code generation for function 'FFTImplementationCallback'
//


// Include files
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <cmath>

// Function Definitions
void FFTImplementationCallback::doHalfLengthBluestein(const coder::array<double,
  1U> &x, coder::array<creal_T, 1U> &y, int nrowsx, int nRows, int nfft, const
  coder::array<creal_T, 1U> &wwc, const coder::array<double, 2U> &costab, const
  coder::array<double, 2U> &sintab, const coder::array<double, 2U> &costabinv,
  const coder::array<double, 2U> &sintabinv)
{
  int hnRows;
  coder::array<creal_T, 1U> ytmp;
  int ju;
  boolean_T tst;
  int ihi;
  int nd2;
  double twid_re;
  int j;
  coder::array<double, 2U> costab1q;
  int k;
  coder::array<double, 2U> b_costab;
  coder::array<double, 2U> b_sintab;
  coder::array<double, 2U> hsintab;
  coder::array<double, 2U> hcostabinv;
  coder::array<double, 2U> hsintabinv;
  int i;
  coder::array<creal_T, 1U> reconVar1;
  int iDelta2;
  coder::array<creal_T, 1U> reconVar2;
  coder::array<int, 2U> wrapIndex;
  int ix;
  int temp_re_tmp;
  double twid_im;
  int loop_ub_tmp;
  coder::array<creal_T, 1U> fy;
  int nRowsD2;
  double temp_re;
  double temp_im;
  coder::array<creal_T, 1U> fv;
  hnRows = nRows / 2;
  ytmp.set_size(hnRows);
  if (hnRows > nrowsx) {
    ytmp.set_size(hnRows);
    for (ju = 0; ju < hnRows; ju++) {
      ytmp[ju].re = 0.0;
      ytmp[ju].im = 0.0;
    }
  }

  if ((x.size(0) & 1) == 0) {
    tst = true;
    ihi = x.size(0);
  } else if (x.size(0) >= nRows) {
    tst = true;
    ihi = nRows;
  } else {
    tst = false;
    ihi = x.size(0) - 1;
  }

  if (ihi >= nRows) {
    ihi = nRows;
  }

  nd2 = nRows << 1;
  twid_re = 6.2831853071795862 / static_cast<double>(nd2);
  j = nd2 / 2 / 2;
  costab1q.set_size(1, (j + 1));
  costab1q[0] = 1.0;
  nd2 = j / 2 - 1;
  for (k = 0; k <= nd2; k++) {
    costab1q[k + 1] = std::cos(twid_re * (static_cast<double>(k) + 1.0));
  }

  ju = nd2 + 2;
  nd2 = j - 1;
  for (k = ju; k <= nd2; k++) {
    costab1q[k] = std::sin(twid_re * static_cast<double>(j - k));
  }

  costab1q[j] = 0.0;
  j = costab1q.size(1) - 1;
  nd2 = (costab1q.size(1) - 1) << 1;
  b_costab.set_size(1, (nd2 + 1));
  b_sintab.set_size(1, (nd2 + 1));
  b_costab[0] = 1.0;
  b_sintab[0] = 0.0;
  for (k = 0; k < j; k++) {
    b_costab[k + 1] = costab1q[k + 1];
    b_sintab[k + 1] = -costab1q[(j - k) - 1];
  }

  ju = costab1q.size(1);
  for (k = ju; k <= nd2; k++) {
    b_costab[k] = -costab1q[nd2 - k];
    b_sintab[k] = -costab1q[k - j];
  }

  nd2 = costab.size(1) / 2;
  costab1q.set_size(1, nd2);
  hsintab.set_size(1, nd2);
  hcostabinv.set_size(1, nd2);
  hsintabinv.set_size(1, nd2);
  for (i = 0; i < nd2; i++) {
    iDelta2 = ((i + 1) << 1) - 2;
    costab1q[i] = costab[iDelta2];
    hsintab[i] = sintab[iDelta2];
    hcostabinv[i] = costabinv[iDelta2];
    hsintabinv[i] = sintabinv[iDelta2];
  }

  reconVar1.set_size(hnRows);
  reconVar2.set_size(hnRows);
  nd2 = 0;
  wrapIndex.set_size(1, hnRows);
  for (i = 0; i < hnRows; i++) {
    reconVar1[i].re = b_sintab[nd2] + 1.0;
    reconVar1[i].im = -b_costab[nd2];
    reconVar2[i].re = 1.0 - b_sintab[nd2];
    reconVar2[i].im = b_costab[nd2];
    nd2 += 2;
    if (i + 1 != 1) {
      wrapIndex[i] = (hnRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }

  nd2 = 0;
  ju = static_cast<int>(static_cast<double>(ihi) / 2.0);
  for (ix = 0; ix < ju; ix++) {
    temp_re_tmp = (hnRows + ix) - 1;
    twid_im = x[nd2 + 1];
    ytmp[ix].re = wwc[temp_re_tmp].re * x[nd2] + wwc[temp_re_tmp].im * twid_im;
    ytmp[ix].im = wwc[temp_re_tmp].re * twid_im - wwc[temp_re_tmp].im * x[nd2];
    nd2 += 2;
  }

  if (!tst) {
    temp_re_tmp = (hnRows + ju) - 1;
    ytmp[ju].re = wwc[temp_re_tmp].re * x[nd2] + wwc[temp_re_tmp].im * 0.0;
    ytmp[ju].im = wwc[temp_re_tmp].re * 0.0 - wwc[temp_re_tmp].im * x[nd2];
    if (ju + 2 <= hnRows) {
      ju = static_cast<int>(static_cast<double>(ihi) / 2.0) + 2;
      for (i = ju; i <= hnRows; i++) {
        ytmp[i - 1].re = 0.0;
        ytmp[i - 1].im = 0.0;
      }
    }
  } else {
    if (ju + 1 <= hnRows) {
      ju = static_cast<int>(static_cast<double>(ihi) / 2.0) + 1;
      for (i = ju; i <= hnRows; i++) {
        ytmp[i - 1].re = 0.0;
        ytmp[i - 1].im = 0.0;
      }
    }
  }

  loop_ub_tmp = static_cast<int>(static_cast<double>(nfft) / 2.0);
  fy.set_size(loop_ub_tmp);
  if (loop_ub_tmp > ytmp.size(0)) {
    fy.set_size(loop_ub_tmp);
    for (ju = 0; ju < loop_ub_tmp; ju++) {
      fy[ju].re = 0.0;
      fy[ju].im = 0.0;
    }
  }

  ihi = ytmp.size(0);
  if (ihi >= loop_ub_tmp) {
    ihi = loop_ub_tmp;
  }

  iDelta2 = loop_ub_tmp - 2;
  nRowsD2 = loop_ub_tmp / 2;
  k = nRowsD2 / 2;
  ix = 0;
  nd2 = 0;
  ju = 0;
  for (i = 0; i <= ihi - 2; i++) {
    fy[nd2] = ytmp[ix];
    j = loop_ub_tmp;
    tst = true;
    while (tst) {
      j >>= 1;
      ju ^= j;
      tst = ((ju & j) == 0);
    }

    nd2 = ju;
    ix++;
  }

  fy[nd2] = ytmp[ix];
  if (loop_ub_tmp > 1) {
    for (i = 0; i <= iDelta2; i += 2) {
      temp_re = fy[i + 1].re;
      temp_im = fy[i + 1].im;
      twid_re = fy[i].re;
      twid_im = fy[i].im;
      fy[i + 1].re = fy[i].re - fy[i + 1].re;
      fy[i + 1].im = fy[i].im - fy[i + 1].im;
      twid_re += temp_re;
      twid_im += temp_im;
      fy[i].re = twid_re;
      fy[i].im = twid_im;
    }
  }

  nd2 = 2;
  iDelta2 = 4;
  ix = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (i = 0; i < ix; i += iDelta2) {
      temp_re_tmp = i + nd2;
      temp_re = fy[temp_re_tmp].re;
      temp_im = fy[temp_re_tmp].im;
      fy[temp_re_tmp].re = fy[i].re - fy[temp_re_tmp].re;
      fy[temp_re_tmp].im = fy[i].im - fy[temp_re_tmp].im;
      fy[i].re = fy[i].re + temp_re;
      fy[i].im = fy[i].im + temp_im;
    }

    ju = 1;
    for (j = k; j < nRowsD2; j += k) {
      twid_re = costab1q[j];
      twid_im = hsintab[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re_tmp = i + nd2;
        temp_re = twid_re * fy[temp_re_tmp].re - twid_im * fy[temp_re_tmp].im;
        temp_im = twid_re * fy[temp_re_tmp].im + twid_im * fy[temp_re_tmp].re;
        fy[temp_re_tmp].re = fy[i].re - temp_re;
        fy[temp_re_tmp].im = fy[i].im - temp_im;
        fy[i].re = fy[i].re + temp_re;
        fy[i].im = fy[i].im + temp_im;
        i += iDelta2;
      }

      ju++;
    }

    k /= 2;
    nd2 = iDelta2;
    iDelta2 += iDelta2;
    ix -= nd2;
  }

  FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (loop_ub_tmp),
    (costab1q), (hsintab), (fv));
  nd2 = fy.size(0);
  for (ju = 0; ju < nd2; ju++) {
    twid_im = fy[ju].re * fv[ju].im + fy[ju].im * fv[ju].re;
    fy[ju].re = fy[ju].re * fv[ju].re - fy[ju].im * fv[ju].im;
    fy[ju].im = twid_im;
  }

  FFTImplementationCallback::r2br_r2dit_trig_impl((fy), (loop_ub_tmp),
    (hcostabinv), (hsintabinv), (fv));
  if (fv.size(0) > 1) {
    twid_re = 1.0 / static_cast<double>(fv.size(0));
    nd2 = fv.size(0);
    for (ju = 0; ju < nd2; ju++) {
      fv[ju].re = twid_re * fv[ju].re;
      fv[ju].im = twid_re * fv[ju].im;
    }
  }

  nd2 = 0;
  ju = wwc.size(0);
  for (k = hnRows; k <= ju; k++) {
    ytmp[nd2].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].im;
    ytmp[nd2].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].re;
    nd2++;
  }

  for (i = 0; i < hnRows; i++) {
    ju = wrapIndex[i];
    twid_re = ytmp[ju - 1].re;
    twid_im = -ytmp[ju - 1].im;
    y[i].re = 0.5 * ((ytmp[i].re * reconVar1[i].re - ytmp[i].im * reconVar1[i].
                      im) + (twid_re * reconVar2[i].re - twid_im * reconVar2[i].
      im));
    y[i].im = 0.5 * ((ytmp[i].re * reconVar1[i].im + ytmp[i].im * reconVar1[i].
                      re) + (twid_re * reconVar2[i].im + twid_im * reconVar2[i].
      re));
    twid_re = ytmp[ju - 1].re;
    twid_im = -ytmp[ju - 1].im;
    ju = hnRows + i;
    y[ju].re = 0.5 * ((ytmp[i].re * reconVar2[i].re - ytmp[i].im * reconVar2[i].
                       im) + (twid_re * reconVar1[i].re - twid_im * reconVar1[i]
      .im));
    y[ju].im = 0.5 * ((ytmp[i].re * reconVar2[i].im + ytmp[i].im * reconVar2[i].
                       re) + (twid_re * reconVar1[i].im + twid_im * reconVar1[i]
      .re));
  }
}

void FFTImplementationCallback::doHalfLengthRadix2(const coder::array<double, 1U>
  &x, coder::array<creal_T, 1U> &y, int unsigned_nRows, const coder::array<
  double, 2U> &costab, const coder::array<double, 2U> &sintab)
{
  int nRows;
  int j;
  int ihi;
  int nRowsD2;
  int k;
  int hszCostab;
  coder::array<double, 2U> hcostab;
  coder::array<double, 2U> hsintab;
  int i;
  coder::array<creal_T, 1U> reconVar1;
  int ju;
  coder::array<creal_T, 1U> reconVar2;
  coder::array<int, 2U> wrapIndex;
  double temp_re;
  double temp_im;
  int loop_ub;
  coder::array<int, 1U> bitrevIndex;
  int iheight;
  boolean_T tst;
  double y_im;
  double y_re;
  double b_y_im;
  double temp2_re;
  double temp2_im;
  nRows = unsigned_nRows / 2;
  j = y.size(0);
  if (j >= nRows) {
    j = nRows;
  }

  ihi = nRows - 2;
  nRowsD2 = nRows / 2;
  k = nRowsD2 / 2;
  hszCostab = costab.size(1) / 2;
  hcostab.set_size(1, hszCostab);
  hsintab.set_size(1, hszCostab);
  for (i = 0; i < hszCostab; i++) {
    ju = ((i + 1) << 1) - 2;
    hcostab[i] = costab[ju];
    hsintab[i] = sintab[ju];
  }

  reconVar1.set_size(nRows);
  reconVar2.set_size(nRows);
  wrapIndex.set_size(1, nRows);
  for (i = 0; i < nRows; i++) {
    temp_re = sintab[i];
    temp_im = costab[i];
    reconVar1[i].re = temp_re + 1.0;
    reconVar1[i].im = -temp_im;
    reconVar2[i].re = 1.0 - temp_re;
    reconVar2[i].im = temp_im;
    if (i + 1 != 1) {
      wrapIndex[i] = (nRows - i) + 1;
    } else {
      wrapIndex[0] = 1;
    }
  }

  ju = 0;
  hszCostab = 1;
  loop_ub = static_cast<int>(static_cast<double>(unsigned_nRows) / 2.0);
  bitrevIndex.set_size(loop_ub);
  for (iheight = 0; iheight < loop_ub; iheight++) {
    bitrevIndex[iheight] = 0;
  }

  for (iheight = 0; iheight <= j - 2; iheight++) {
    bitrevIndex[iheight] = hszCostab;
    hszCostab = loop_ub;
    tst = true;
    while (tst) {
      hszCostab >>= 1;
      ju ^= hszCostab;
      tst = ((ju & hszCostab) == 0);
    }

    hszCostab = ju + 1;
  }

  bitrevIndex[j - 1] = hszCostab;
  if ((x.size(0) & 1) == 0) {
    tst = true;
    j = x.size(0);
  } else if (x.size(0) >= unsigned_nRows) {
    tst = true;
    j = unsigned_nRows;
  } else {
    tst = false;
    j = x.size(0) - 1;
  }

  hszCostab = 0;
  if (j >= unsigned_nRows) {
    j = unsigned_nRows;
  }

  temp_re = static_cast<double>(j) / 2.0;
  iheight = static_cast<int>(temp_re);
  for (i = 0; i < iheight; i++) {
    y[bitrevIndex[i] - 1].re = x[hszCostab];
    y[bitrevIndex[i] - 1].im = x[hszCostab + 1];
    hszCostab += 2;
  }

  if (!tst) {
    iheight = bitrevIndex[static_cast<int>(temp_re)] - 1;
    y[iheight].re = x[hszCostab];
    y[iheight].im = 0.0;
  }

  if (nRows > 1) {
    for (i = 0; i <= ihi; i += 2) {
      temp_re = y[i + 1].re;
      temp_im = y[i + 1].im;
      y[i + 1].re = y[i].re - y[i + 1].re;
      y[i + 1].im = y[i].im - y[i + 1].im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }
  }

  hszCostab = 2;
  ju = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int temp_re_tmp;
    for (i = 0; i < iheight; i += ju) {
      temp_re_tmp = i + hszCostab;
      temp_re = y[temp_re_tmp].re;
      temp_im = y[temp_re_tmp].im;
      y[temp_re_tmp].re = y[i].re - temp_re;
      y[temp_re_tmp].im = y[i].im - temp_im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }

    nRows = 1;
    for (j = k; j < nRowsD2; j += k) {
      temp2_re = hcostab[j];
      temp2_im = hsintab[j];
      i = nRows;
      ihi = nRows + iheight;
      while (i < ihi) {
        temp_re_tmp = i + hszCostab;
        temp_re = temp2_re * y[temp_re_tmp].re - temp2_im * y[temp_re_tmp].im;
        temp_im = temp2_re * y[temp_re_tmp].im + temp2_im * y[temp_re_tmp].re;
        y[temp_re_tmp].re = y[i].re - temp_re;
        y[temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += ju;
      }

      nRows++;
    }

    k /= 2;
    hszCostab = ju;
    ju += ju;
    iheight -= hszCostab;
  }

  hszCostab = loop_ub / 2;
  temp_re = y[0].re;
  temp_im = y[0].im;
  y_im = y[0].re * reconVar1[0].im + y[0].im * reconVar1[0].re;
  y_re = y[0].re;
  b_y_im = -y[0].im;
  y[0].re = 0.5 * ((y[0].re * reconVar1[0].re - y[0].im * reconVar1[0].im) +
                   (y_re * reconVar2[0].re - b_y_im * reconVar2[0].im));
  y[0].im = 0.5 * (y_im + (y_re * reconVar2[0].im + b_y_im * reconVar2[0].re));
  y[loop_ub].re = 0.5 * ((temp_re * reconVar2[0].re - temp_im * reconVar2[0].im)
    + (temp_re * reconVar1[0].re - -temp_im * reconVar1[0].im));
  y[loop_ub].im = 0.5 * ((temp_re * reconVar2[0].im + temp_im * reconVar2[0].re)
    + (temp_re * reconVar1[0].im + -temp_im * reconVar1[0].re));
  for (i = 2; i <= hszCostab; i++) {
    temp_re = y[i - 1].re;
    temp_im = y[i - 1].im;
    iheight = wrapIndex[i - 1];
    temp2_re = y[iheight - 1].re;
    temp2_im = y[iheight - 1].im;
    y_im = y[i - 1].re * reconVar1[i - 1].im + y[i - 1].im * reconVar1[i - 1].re;
    y_re = y[iheight - 1].re;
    b_y_im = -y[iheight - 1].im;
    y[i - 1].re = 0.5 * ((y[i - 1].re * reconVar1[i - 1].re - y[i - 1].im *
                          reconVar1[i - 1].im) + (y_re * reconVar2[i - 1].re -
      b_y_im * reconVar2[i - 1].im));
    y[i - 1].im = 0.5 * (y_im + (y_re * reconVar2[i - 1].im + b_y_im *
      reconVar2[i - 1].re));
    ju = (loop_ub + i) - 1;
    y[ju].re = 0.5 * ((temp_re * reconVar2[i - 1].re - temp_im * reconVar2[i - 1]
                       .im) + (temp2_re * reconVar1[i - 1].re - -temp2_im *
      reconVar1[i - 1].im));
    y[ju].im = 0.5 * ((temp_re * reconVar2[i - 1].im + temp_im * reconVar2[i - 1]
                       .re) + (temp2_re * reconVar1[i - 1].im + -temp2_im *
      reconVar1[i - 1].re));
    y[iheight - 1].re = 0.5 * ((temp2_re * reconVar1[iheight - 1].re - temp2_im *
      reconVar1[iheight - 1].im) + (temp_re * reconVar2[iheight - 1].re -
      -temp_im * reconVar2[iheight - 1].im));
    y[iheight - 1].im = 0.5 * ((temp2_re * reconVar1[iheight - 1].im + temp2_im *
      reconVar1[iheight - 1].re) + (temp_re * reconVar2[iheight - 1].im +
      -temp_im * reconVar2[iheight - 1].re));
    ju = (iheight + loop_ub) - 1;
    y[ju].re = 0.5 * ((temp2_re * reconVar2[iheight - 1].re - temp2_im *
                       reconVar2[iheight - 1].im) + (temp_re * reconVar1[iheight
      - 1].re - -temp_im * reconVar1[iheight - 1].im));
    y[ju].im = 0.5 * ((temp2_re * reconVar2[iheight - 1].im + temp2_im *
                       reconVar2[iheight - 1].re) + (temp_re * reconVar1[iheight
      - 1].im + -temp_im * reconVar1[iheight - 1].re));
  }

  if (hszCostab != 0) {
    temp_re = y[hszCostab].re;
    temp_im = y[hszCostab].im;
    y_im = y[hszCostab].re * reconVar1[hszCostab].im + y[hszCostab].im *
      reconVar1[hszCostab].re;
    y_re = y[hszCostab].re;
    b_y_im = -y[hszCostab].im;
    y[hszCostab].re = 0.5 * ((y[hszCostab].re * reconVar1[hszCostab].re -
      y[hszCostab].im * reconVar1[hszCostab].im) + (y_re * reconVar2[hszCostab].
      re - b_y_im * reconVar2[hszCostab].im));
    y[hszCostab].im = 0.5 * (y_im + (y_re * reconVar2[hszCostab].im + b_y_im *
      reconVar2[hszCostab].re));
    iheight = loop_ub + hszCostab;
    y[iheight].re = 0.5 * ((temp_re * reconVar2[hszCostab].re - temp_im *
      reconVar2[hszCostab].im) + (temp_re * reconVar1[hszCostab].re - -temp_im *
      reconVar1[hszCostab].im));
    y[iheight].im = 0.5 * ((temp_re * reconVar2[hszCostab].im + temp_im *
      reconVar2[hszCostab].re) + (temp_re * reconVar1[hszCostab].im + -temp_im *
      reconVar1[hszCostab].re));
  }
}

void FFTImplementationCallback::dobluesteinfft(const coder::array<double, 1U> &x,
  int n2blue, int nfft, const coder::array<double, 2U> &costab, const coder::
  array<double, 2U> &sintab, const coder::array<double, 2U> &sintabinv, coder::
  array<creal_T, 1U> &y)
{
  int nInt2m1;
  coder::array<creal_T, 1U> wwc;
  int idx;
  int rt;
  int nInt2;
  int k;
  double nt_re;
  coder::array<creal_T, 1U> fv;
  coder::array<creal_T, 1U> b_fv;
  if ((nfft != 1) && ((nfft & 1) == 0)) {
    int nRows;
    nRows = nfft / 2;
    nInt2m1 = (nRows + nRows) - 1;
    wwc.set_size(nInt2m1);
    idx = nRows;
    rt = 0;
    wwc[nRows - 1].re = 1.0;
    wwc[nRows - 1].im = 0.0;
    nInt2 = nRows << 1;
    for (k = 0; k <= nRows - 2; k++) {
      int b_y;
      double nt_im;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }

      nt_im = -3.1415926535897931 * static_cast<double>(rt) / static_cast<double>
        (nRows);
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
    rt = nInt2m1 - 1;
    for (k = rt; k >= nRows; k--) {
      wwc[k] = wwc[idx];
      idx++;
    }
  } else {
    nInt2m1 = (nfft + nfft) - 1;
    wwc.set_size(nInt2m1);
    idx = nfft;
    rt = 0;
    wwc[nfft - 1].re = 1.0;
    wwc[nfft - 1].im = 0.0;
    nInt2 = nfft << 1;
    for (k = 0; k <= nfft - 2; k++) {
      int b_y;
      double nt_im;
      b_y = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= b_y) {
        rt += b_y - nInt2;
      } else {
        rt += b_y;
      }

      nt_im = -3.1415926535897931 * static_cast<double>(rt) / static_cast<double>
        (nfft);
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
    rt = nInt2m1 - 1;
    for (k = rt; k >= nfft; k--) {
      wwc[k] = wwc[idx];
      idx++;
    }
  }

  y.set_size(nfft);
  if (nfft > x.size(0)) {
    y.set_size(nfft);
    for (rt = 0; rt < nfft; rt++) {
      y[rt].re = 0.0;
      y[rt].im = 0.0;
    }
  }

  if ((n2blue != 1) && ((nfft & 1) == 0)) {
    FFTImplementationCallback::doHalfLengthBluestein((x), (y), (x.size(0)),
      (nfft), (n2blue), (wwc), (costab), (sintab), (costab), (sintabinv));
  } else {
    nInt2m1 = x.size(0);
    if (nfft < nInt2m1) {
      nInt2m1 = nfft;
    }

    rt = 0;
    for (k = 0; k < nInt2m1; k++) {
      nInt2 = (nfft + k) - 1;
      y[k].re = wwc[nInt2].re * x[rt];
      y[k].im = wwc[nInt2].im * -x[rt];
      rt++;
    }

    rt = nInt2m1 + 1;
    for (k = rt; k <= nfft; k++) {
      y[k - 1].re = 0.0;
      y[k - 1].im = 0.0;
    }

    FFTImplementationCallback::r2br_r2dit_trig_impl((y), (n2blue), (costab),
      (sintab), (fv));
    FFTImplementationCallback::r2br_r2dit_trig_impl((wwc), (n2blue), (costab),
      (sintab), (b_fv));
    b_fv.set_size(fv.size(0));
    nInt2m1 = fv.size(0);
    for (rt = 0; rt < nInt2m1; rt++) {
      nt_re = fv[rt].re * b_fv[rt].im + fv[rt].im * b_fv[rt].re;
      b_fv[rt].re = fv[rt].re * b_fv[rt].re - fv[rt].im * b_fv[rt].im;
      b_fv[rt].im = nt_re;
    }

    FFTImplementationCallback::r2br_r2dit_trig_impl((b_fv), (n2blue), (costab),
      (sintabinv), (fv));
    if (fv.size(0) > 1) {
      nt_re = 1.0 / static_cast<double>(fv.size(0));
      nInt2m1 = fv.size(0);
      for (rt = 0; rt < nInt2m1; rt++) {
        fv[rt].re = nt_re * fv[rt].re;
        fv[rt].im = nt_re * fv[rt].im;
      }
    }

    idx = 0;
    rt = wwc.size(0);
    for (k = nfft; k <= rt; k++) {
      y[idx].re = wwc[k - 1].re * fv[k - 1].re + wwc[k - 1].im * fv[k - 1].im;
      y[idx].im = wwc[k - 1].re * fv[k - 1].im - wwc[k - 1].im * fv[k - 1].re;
      idx++;
    }
  }
}

void FFTImplementationCallback::get_algo_sizes(int nfft, boolean_T useRadix2,
  int *n2blue, int *nRows)
{
  *n2blue = 1;
  if (useRadix2) {
    *nRows = nfft;
  } else {
    if (nfft > 0) {
      int n;
      int pmax;
      n = (nfft + nfft) - 1;
      pmax = 31;
      if (n <= 1) {
        pmax = 0;
      } else {
        int pmin;
        boolean_T exitg1;
        pmin = 0;
        exitg1 = false;
        while ((!exitg1) && (pmax - pmin > 1)) {
          int k;
          int pow2p;
          k = (pmin + pmax) >> 1;
          pow2p = 1 << k;
          if (pow2p == n) {
            pmax = k;
            exitg1 = true;
          } else if (pow2p > n) {
            pmax = k;
          } else {
            pmin = k;
          }
        }
      }

      *n2blue = 1 << pmax;
    }

    *nRows = *n2blue;
  }
}

void FFTImplementationCallback::r2br_r2dit_trig_impl(const coder::array<creal_T,
  1U> &x, int unsigned_nRows, const coder::array<double, 2U> &costab, const
  coder::array<double, 2U> &sintab, coder::array<creal_T, 1U> &y)
{
  int iDelta2;
  int iy;
  int iheight;
  int nRowsD2;
  int k;
  int ix;
  int ju;
  int i;
  double temp_re;
  double temp_im;
  double twid_re;
  double twid_im;
  y.set_size(unsigned_nRows);
  if (unsigned_nRows > x.size(0)) {
    y.set_size(unsigned_nRows);
    for (iy = 0; iy < unsigned_nRows; iy++) {
      y[iy].re = 0.0;
      y[iy].im = 0.0;
    }
  }

  iDelta2 = x.size(0);
  if (iDelta2 >= unsigned_nRows) {
    iDelta2 = unsigned_nRows;
  }

  iheight = unsigned_nRows - 2;
  nRowsD2 = unsigned_nRows / 2;
  k = nRowsD2 / 2;
  ix = 0;
  iy = 0;
  ju = 0;
  for (i = 0; i <= iDelta2 - 2; i++) {
    boolean_T tst;
    y[iy] = x[ix];
    iy = unsigned_nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= iheight; i += 2) {
      temp_re = y[i + 1].re;
      temp_im = y[i + 1].im;
      twid_re = y[i].re;
      twid_im = y[i].im;
      y[i + 1].re = y[i].re - y[i + 1].re;
      y[i + 1].im = y[i].im - y[i + 1].im;
      twid_re += temp_re;
      twid_im += temp_im;
      y[i].re = twid_re;
      y[i].im = twid_im;
    }
  }

  iy = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int temp_re_tmp;
    for (i = 0; i < iheight; i += iDelta2) {
      temp_re_tmp = i + iy;
      temp_re = y[temp_re_tmp].re;
      temp_im = y[temp_re_tmp].im;
      y[temp_re_tmp].re = y[i].re - y[temp_re_tmp].re;
      y[temp_re_tmp].im = y[i].im - y[temp_re_tmp].im;
      y[i].re = y[i].re + temp_re;
      y[i].im = y[i].im + temp_im;
    }

    ix = 1;
    for (ju = k; ju < nRowsD2; ju += k) {
      int ihi;
      twid_re = costab[ju];
      twid_im = sintab[ju];
      i = ix;
      ihi = ix + iheight;
      while (i < ihi) {
        temp_re_tmp = i + iy;
        temp_re = twid_re * y[temp_re_tmp].re - twid_im * y[temp_re_tmp].im;
        temp_im = twid_re * y[temp_re_tmp].im + twid_im * y[temp_re_tmp].re;
        y[temp_re_tmp].re = y[i].re - temp_re;
        y[temp_re_tmp].im = y[i].im - temp_im;
        y[i].re = y[i].re + temp_re;
        y[i].im = y[i].im + temp_im;
        i += iDelta2;
      }

      ix++;
    }

    k /= 2;
    iy = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iy;
  }
}

// End of code generation (FFTImplementationCallback.cpp)
