//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  composite_cpp.cpp
//
//  Code generation for function 'composite_cpp'
//


// Include files
#include <composite_cpp.h>
#include <FFTImplementationCallback.h>
#include <composite_cpp_data.h>
#include <composite_cpp_initialize.h>
#include <composite_cpp_rtwutil.h>
#include <fft.h>
#include <gencoswin.h>
#include <hann.h>
#include <minOrMax.h>
#include <pesq_original_cpp.h>
#include <rt_nonfinite.h>
#include <sort.h>
#include <cmath>
#include <cstring>
#include <math.h>

// Function Declarations
static void llr(const coder::array<double, 1U> &clean_speech, const coder::array<
                double, 1U> &processed_speech, double sample_rate, coder::array<
                double, 2U> &distortion);
static void lpcoeff(const coder::array<double, 1U> &speech_frame, double
                    model_order, double acorr_data[], int acorr_size[2], double
                    refcoeff_data[], int refcoeff_size[2], double lpparams_data[],
                    int lpparams_size[2]);
static void wss(const coder::array<double, 1U> &clean_speech, const coder::array<
                double, 1U> &processed_speech, double sample_rate, coder::array<
                double, 2U> &distortion);

// Function Definitions
static void llr(const coder::array<double, 1U> &clean_speech, const coder::array<
                double, 1U> &processed_speech, double sample_rate, coder::array<
                double, 2U> &distortion)
{
  double winlength;
  double skiprate;
  int P;
  double num_frames;
  double start;
  coder::array<double, 2U> r;
  int nx;
  coder::array<double, 1U> window;
  int i;
  int k;
  coder::array<double, 1U> b_clean_speech;
  double R_clean_data[17];
  int R_clean_size[2];
  double Ref_clean_data[16];
  int Ref_clean_size[2];
  double A_clean_data[17];
  int A_clean_size[2];
  double R_processed_data[17];
  int R_processed_size[2];
  double A_processed_data[17];
  int A_processed_size[2];
  double b_tmp_data[289];
  double R_processed;

  // -----------------------------------------------
  //  ----------------------------------------------------------------------
  //  Check the length of the clean and processed speech.  Must be the same.
  //  ----------------------------------------------------------------------
  //  if (clean_length ~= processed_length)
  //    disp('Error: Both Speech Files must be same length.');
  //    return
  //  end
  //  ----------------------------------------------------------------------
  //  Global Variables
  //  ----------------------------------------------------------------------
  winlength = rt_roundd_snf(30.0 * sample_rate / 1000.0);

  //   window length in samples
  skiprate = std::floor(winlength / 4.0);

  //  window skip in samples
  if (sample_rate < 10000.0) {
    P = 10;

    //  LPC Analysis Order
  } else {
    P = 16;

    //  this could vary depending on sampling frequency.
  }

  //  ----------------------------------------------------------------------
  //  For each frame of input speech, calculate the Log Likelihood Ratio
  //  ----------------------------------------------------------------------
  num_frames = static_cast<double>(clean_speech.size(0)) / skiprate - winlength /
    skiprate;

  //  number of frames
  start = 1.0;

  //  starting sample
  if (winlength < 1.0) {
    r.set_size(1, 0);
  } else if (rtIsInf(winlength) && (1.0 == winlength)) {
    r.set_size(1, 1);
    r[0] = rtNaN;
  } else {
    nx = static_cast<int>(winlength - 1.0);
    r.set_size(1, (nx + 1));
    for (i = 0; i <= nx; i++) {
      r[i] = static_cast<double>(i) + 1.0;
    }
  }

  window.set_size(r.size(1));
  nx = r.size(1);
  for (i = 0; i < nx; i++) {
    window[i] = 6.2831853071795862 * r[i] / (winlength + 1.0);
  }

  nx = window.size(0);
  for (k = 0; k < nx; k++) {
    window[k] = std::cos(window[k]);
  }

  nx = window.size(0);
  for (i = 0; i < nx; i++) {
    window[i] = 0.5 * (1.0 - window[i]);
  }

  nx = static_cast<int>(std::floor(num_frames));
  distortion.set_size(1, nx);
  for (i = 0; i < nx; i++) {
    distortion[i] = 0.0;
  }

  i = static_cast<int>(num_frames);
  for (int frame_count = 0; frame_count < i; frame_count++) {
    int i1;
    int i2;
    int boffset;
    int b_tmp_size_idx_0;
    int b_tmp_size_idx_1;
    int j;

    //  ----------------------------------------------------------
    //  (1) Get the Frames for the test and reference speech.
    //      Multiply by Hanning Window.
    //  ----------------------------------------------------------
    num_frames = (start + winlength) - 1.0;
    if (start > num_frames) {
      i1 = 0;
      i2 = 0;
      boffset = 0;
      b_tmp_size_idx_0 = 0;
    } else {
      i1 = static_cast<int>(start) - 1;
      i2 = static_cast<int>(num_frames);
      boffset = i1;
      b_tmp_size_idx_0 = i2;
    }

    //  ----------------------------------------------------------
    //  (2) Get the autocorrelation lags and LPC parameters used
    //      to compute the LLR measure.
    //  ----------------------------------------------------------
    nx = i2 - i1;
    b_clean_speech.set_size(nx);
    for (i2 = 0; i2 < nx; i2++) {
      b_clean_speech[i2] = clean_speech[i1 + i2] * window[i2];
    }

    lpcoeff(b_clean_speech, static_cast<double>(P), R_clean_data, R_clean_size,
            Ref_clean_data, Ref_clean_size, A_clean_data, A_clean_size);
    nx = b_tmp_size_idx_0 - boffset;
    b_clean_speech.set_size(nx);
    for (i1 = 0; i1 < nx; i1++) {
      b_clean_speech[i1] = processed_speech[boffset + i1] * window[i1];
    }

    lpcoeff(b_clean_speech, static_cast<double>(P), R_processed_data,
            R_processed_size, Ref_clean_data, Ref_clean_size, A_processed_data,
            A_processed_size);

    //  ----------------------------------------------------------
    //  (3) Compute the LLR measure
    //  ----------------------------------------------------------
    b_tmp_size_idx_0 = static_cast<signed char>(R_clean_size[1]);
    b_tmp_size_idx_1 = static_cast<signed char>(R_clean_size[1]);
    nx = 0;
    i1 = R_clean_size[1] - 1;
    for (j = 0; j <= i1; j++) {
      k = j;
      i2 = R_clean_size[1] - 1;
      for (boffset = 0; boffset <= i2; boffset++) {
        if (boffset < j) {
          b_tmp_data[nx] = R_clean_data[k];
          k--;
        } else {
          b_tmp_data[nx] = R_clean_data[k];
          k++;
        }

        nx++;
      }
    }

    nx = A_processed_size[1];
    for (j = 0; j < b_tmp_size_idx_1; j++) {
      boffset = j * b_tmp_size_idx_0;
      R_clean_data[j] = 0.0;
      for (k = 0; k < nx; k++) {
        R_clean_data[j] += A_processed_data[k] * b_tmp_data[boffset + k];
      }
    }

    nx = A_clean_size[1];
    for (j = 0; j < b_tmp_size_idx_1; j++) {
      boffset = j * b_tmp_size_idx_0;
      R_processed_data[j] = 0.0;
      for (k = 0; k < nx; k++) {
        R_processed_data[j] += A_clean_data[k] * b_tmp_data[boffset + k];
      }
    }

    num_frames = 0.0;
    for (i1 = 0; i1 < b_tmp_size_idx_1; i1++) {
      num_frames += R_clean_data[i1] * A_processed_data[i1];
    }

    R_processed = 0.0;
    for (i1 = 0; i1 < b_tmp_size_idx_1; i1++) {
      R_processed += R_processed_data[i1] * A_clean_data[i1];
    }

    distortion[frame_count] = std::log(num_frames / R_processed);
    start += skiprate;
  }
}

static void lpcoeff(const coder::array<double, 1U> &speech_frame, double
                    model_order, double acorr_data[], int acorr_size[2], double
                    refcoeff_data[], int refcoeff_size[2], double lpparams_data[],
                    int lpparams_size[2])
{
  int winlength;
  int loop_ub;
  int k;
  int vlen;
  int b_loop_ub;
  double a_data[16];
  double a_past_data[16];
  coder::array<double, 1U> x;
  double E_data[17];
  int i;
  double y;
  int b_k;
  coder::array<double, 2U> b_x;

  // ---------------------------------------------
  //  ----------------------------------------------------------
  //  (1) Compute Autocorrelation Lags
  //  ----------------------------------------------------------
  if (speech_frame.size(0) < 1) {
    winlength = 1;
  } else {
    winlength = speech_frame.size(0);
  }

  acorr_size[0] = 1;
  loop_ub = static_cast<int>(model_order + 1.0);
  acorr_size[1] = loop_ub;
  for (k = 0; k < loop_ub; k++) {
    acorr_data[k] = 0.0;
    vlen = static_cast<int>(static_cast<double>(winlength) - (static_cast<double>
      (k) + 1.0));
    if (1 > vlen + 1) {
      b_loop_ub = -1;
    } else {
      b_loop_ub = vlen;
    }

    if (k + 1U > static_cast<unsigned int>(winlength)) {
      vlen = 0;
    } else {
      vlen = k;
    }

    x.set_size((b_loop_ub + 1));
    for (i = 0; i <= b_loop_ub; i++) {
      x[i] = speech_frame[i] * speech_frame[vlen + i];
    }

    vlen = x.size(0);
    if (x.size(0) == 0) {
      y = 0.0;
    } else {
      y = x[0];
      for (b_k = 2; b_k <= vlen; b_k++) {
        y += x[b_k - 1];
      }
    }

    acorr_data[k] = y;
  }

  //  ----------------------------------------------------------
  //  (2) Levinson-Durbin
  //  ----------------------------------------------------------
  b_loop_ub = static_cast<int>(model_order);
  for (vlen = 0; vlen < b_loop_ub; vlen++) {
    a_data[vlen] = 1.0;
  }

  if (0 <= b_loop_ub - 1) {
    std::memset(&a_past_data[0], 0, b_loop_ub * sizeof(double));
  }

  if (0 <= loop_ub - 1) {
    std::memset(&E_data[0], 0, loop_ub * sizeof(double));
  }

  refcoeff_size[0] = 1;
  refcoeff_size[1] = b_loop_ub;
  E_data[0] = acorr_data[0];
  for (winlength = 0; winlength < b_loop_ub; winlength++) {
    if (1 > winlength) {
      loop_ub = 0;
    } else {
      loop_ub = winlength;
    }

    if (0 <= loop_ub - 1) {
      std::memcpy(&a_past_data[0], &a_data[0], loop_ub * sizeof(double));
    }

    if (1 > winlength) {
      loop_ub = 0;
    } else {
      loop_ub = winlength;
    }

    if (2U > winlength + 1U) {
      vlen = 0;
      i = 1;
    } else {
      vlen = winlength;
      i = -1;
    }

    b_x.set_size(1, loop_ub);
    for (b_k = 0; b_k < loop_ub; b_k++) {
      b_x[b_k] = a_past_data[b_k] * acorr_data[vlen + i * b_k];
    }

    vlen = b_x.size(1);
    if (b_x.size(1) == 0) {
      y = 0.0;
    } else {
      y = b_x[0];
      for (k = 2; k <= vlen; k++) {
        y += b_x[k - 1];
      }
    }

    y = (acorr_data[winlength + 1] - y) / E_data[winlength];
    refcoeff_data[winlength] = y;
    a_data[winlength] = y;
    if (1 > winlength) {
      loop_ub = 0;
      vlen = 1;
      i = 1;
    } else {
      loop_ub = winlength;
      vlen = winlength;
      i = -1;
    }

    for (b_k = 0; b_k < loop_ub; b_k++) {
      a_data[b_k] = a_past_data[b_k] - y * a_past_data[(vlen + i * b_k) - 1];
    }

    E_data[winlength + 1] = (1.0 - refcoeff_data[winlength] *
      refcoeff_data[winlength]) * E_data[winlength];
  }

  lpparams_size[0] = 1;
  lpparams_size[1] = b_loop_ub + 1;
  lpparams_data[0] = 1.0;
  for (vlen = 0; vlen < b_loop_ub; vlen++) {
    lpparams_data[vlen + 1] = -a_data[vlen];
  }
}

static void wss(const coder::array<double, 1U> &clean_speech, const coder::array<
                double, 1U> &processed_speech, double sample_rate, coder::array<
                double, 2U> &distortion)
{
  double winlength;
  double skiprate;
  double max_freq;
  double absn;
  double dBMax_processed;
  int nx;
  double n_fft;
  double n_fftby2;
  double cent_freq[25];
  double bandwidth[25];
  int i;
  coder::array<double, 2U> crit_filter;
  int b_i;
  coder::array<double, 2U> j;
  double x;
  coder::array<double, 1U> window;
  coder::array<double, 2U> y;
  int k;
  int loop_ub;
  int b_loop_ub;
  coder::array<double, 1U> b_y;
  coder::array<creal_T, 1U> b_x;
  coder::array<double, 1U> clean_spec;
  coder::array<double, 1U> processed_spec;
  double clean_slope[24];
  double processed_slope[24];
  double clean_loc_peak[24];

  // =================================================================
  //  ----------------------------------------------------------------------
  //  Check the length of the clean and processed speech.  Must be the same.
  //  ----------------------------------------------------------------------
  //  if (clean_length ~= processed_length)
  //    disp('Error: Files  musthave same length.');
  //    return
  //  end
  //  ----------------------------------------------------------------------
  //  Global Variables
  //  ----------------------------------------------------------------------
  winlength = rt_roundd_snf(30.0 * sample_rate / 1000.0);

  // 240;		   % window length in samples
  skiprate = std::floor(winlength / 4.0);

  //  window skip in samples
  max_freq = sample_rate / 2.0;

  //  maximum bandwidth
  //  number of critical bands
  //  defaults to 10th order LP spectrum
  absn = std::abs(2.0 * winlength);
  if ((!rtIsInf(absn)) && (!rtIsNaN(absn))) {
    dBMax_processed = frexp(absn, &nx);
    absn = nx;
    if (dBMax_processed == 0.5) {
      absn = static_cast<double>(nx) - 1.0;
    }
  }

  n_fft = rt_powd_snf(2.0, absn);
  n_fftby2 = n_fft / 2.0;

  //  FFT size/2
  //  value suggested by Klatt, pg 1280
  //  value suggested by Klatt, pg 1280		
  //  ----------------------------------------------------------------------
  //  Critical Band Filter Definitions (Center Frequency and Bandwidths in Hz)
  //  ----------------------------------------------------------------------
  cent_freq[0] = 50.0;
  bandwidth[0] = 70.0;
  cent_freq[1] = 120.0;
  bandwidth[1] = 70.0;
  cent_freq[2] = 190.0;
  bandwidth[2] = 70.0;
  cent_freq[3] = 260.0;
  bandwidth[3] = 70.0;
  cent_freq[4] = 330.0;
  bandwidth[4] = 70.0;
  cent_freq[5] = 400.0;
  bandwidth[5] = 70.0;
  cent_freq[6] = 470.0;
  bandwidth[6] = 70.0;
  cent_freq[7] = 540.0;
  bandwidth[7] = 77.3724;
  cent_freq[8] = 617.372;
  bandwidth[8] = 86.0056;
  cent_freq[9] = 703.378;
  bandwidth[9] = 95.3398;
  cent_freq[10] = 798.717;
  bandwidth[10] = 105.411;
  cent_freq[11] = 904.128;
  bandwidth[11] = 116.256;
  cent_freq[12] = 1020.38;
  bandwidth[12] = 127.914;
  cent_freq[13] = 1148.3;
  bandwidth[13] = 140.423;
  cent_freq[14] = 1288.72;
  bandwidth[14] = 153.823;
  cent_freq[15] = 1442.54;
  bandwidth[15] = 168.154;
  cent_freq[16] = 1610.7;
  bandwidth[16] = 183.457;
  cent_freq[17] = 1794.16;
  bandwidth[17] = 199.776;
  cent_freq[18] = 1993.93;
  bandwidth[18] = 217.153;
  cent_freq[19] = 2211.08;
  bandwidth[19] = 235.631;
  cent_freq[20] = 2446.71;
  bandwidth[20] = 255.255;
  cent_freq[21] = 2701.97;
  bandwidth[21] = 276.072;
  cent_freq[22] = 2978.04;
  bandwidth[22] = 298.126;
  cent_freq[23] = 3276.17;
  bandwidth[23] = 321.465;
  cent_freq[24] = 3597.63;
  bandwidth[24] = 346.136;

  //  minimum critical bandwidth
  //  ----------------------------------------------------------------------
  //  Set up the critical band filters.  Note here that Gaussianly shaped
  //  filters are used.  Also, the sum of the filter weights are equivalent
  //  for each critical band filter.  Filter less than -30 dB and set to
  //  zero.
  //  ----------------------------------------------------------------------
  //  -30 dB point of filter
  i = static_cast<int>(n_fftby2);
  crit_filter.set_size(25, i);
  nx = 25 * i;
  for (i = 0; i < nx; i++) {
    crit_filter[i] = 0.0;
  }

  for (b_i = 0; b_i < 25; b_i++) {
    absn = bandwidth[b_i] / max_freq * n_fftby2;
    dBMax_processed = std::log(bandwidth[b_i]);
    if (n_fftby2 - 1.0 < 0.0) {
      j.set_size(1, 0);
    } else if (rtIsInf(n_fftby2 - 1.0) && (0.0 == n_fftby2 - 1.0)) {
      j.set_size(1, 1);
      j[0] = rtNaN;
    } else {
      nx = static_cast<int>(std::floor(n_fftby2 - 1.0));
      j.set_size(1, (nx + 1));
      for (i = 0; i <= nx; i++) {
        j[i] = i;
      }
    }

    x = std::floor(cent_freq[b_i] / max_freq * n_fftby2);
    i = j.size(0) * j.size(1);
    j.set_size(1, j.size(1));
    nx = i - 1;
    for (i = 0; i <= nx; i++) {
      j[i] = (j[i] - x) / absn;
    }

    y.set_size(1, j.size(1));
    nx = j.size(1);
    for (k = 0; k < nx; k++) {
      y[k] = rt_powd_snf(j[k], 2.0);
    }

    i = y.size(0) * y.size(1);
    y.set_size(1, y.size(1));
    nx = i - 1;
    for (i = 0; i <= nx; i++) {
      y[i] = -11.0 * y[i] + (4.2484952420493594 - dBMax_processed);
    }

    nx = y.size(1);
    for (k = 0; k < nx; k++) {
      y[k] = std::exp(y[k]);
    }

    nx = y.size(1);
    for (i = 0; i < nx; i++) {
      crit_filter[b_i + 25 * i] = y[i];
    }

    nx = crit_filter.size(1) - 1;
    j.set_size(1, crit_filter.size(1));
    for (i = 0; i <= nx; i++) {
      x = crit_filter[b_i + 25 * i];
      j[i] = x * static_cast<double>(x > 0.0014836595188319822);
    }

    nx = j.size(1);
    for (i = 0; i < nx; i++) {
      crit_filter[b_i + 25 * i] = j[i];
    }
  }

  //  ----------------------------------------------------------------------
  //  For each frame of input speech, calculate the Weighted Spectral
  //  Slope Measure
  //  ----------------------------------------------------------------------
  absn = static_cast<double>(clean_speech.size(0)) / skiprate - winlength /
    skiprate;

  //  number of frames
  max_freq = 1.0;

  //  starting sample
  if (winlength < 1.0) {
    j.set_size(1, 0);
  } else if (rtIsInf(winlength) && (1.0 == winlength)) {
    j.set_size(1, 1);
    j[0] = rtNaN;
  } else {
    nx = static_cast<int>(winlength - 1.0);
    j.set_size(1, (nx + 1));
    for (i = 0; i <= nx; i++) {
      j[i] = static_cast<double>(i) + 1.0;
    }
  }

  window.set_size(j.size(1));
  nx = j.size(1);
  for (i = 0; i < nx; i++) {
    window[i] = 6.2831853071795862 * j[i] / (winlength + 1.0);
  }

  nx = window.size(0);
  for (k = 0; k < nx; k++) {
    window[k] = std::cos(window[k]);
  }

  nx = window.size(0);
  for (i = 0; i < nx; i++) {
    window[i] = 0.5 * (1.0 - window[i]);
  }

  nx = static_cast<int>(std::floor(absn));
  distortion.set_size(1, nx);
  for (i = 0; i < nx; i++) {
    distortion[i] = 0.0;
  }

  i = static_cast<int>(absn);
  if (0 <= i - 1) {
    if (1.0 > n_fftby2) {
      loop_ub = 0;
      b_loop_ub = 0;
    } else {
      loop_ub = static_cast<int>(n_fftby2);
      b_loop_ub = static_cast<int>(n_fftby2);
    }
  }

  for (int frame_count = 0; frame_count < i; frame_count++) {
    int i1;
    int i2;

    //  ----------------------------------------------------------
    //  (1) Get the Frames for the test and reference speech.
    //      Multiply by Hanning Window.
    //  ----------------------------------------------------------
    x = (max_freq + winlength) - 1.0;
    if (max_freq > x) {
      i1 = 0;
      k = 0;
      b_i = 0;
      i2 = 0;
    } else {
      i1 = static_cast<int>(max_freq) - 1;
      k = static_cast<int>(x);
      b_i = i1;
      i2 = k;
    }

    //  ----------------------------------------------------------
    //  (2) Compute the Power Spectrum of Clean and Processed
    //  ----------------------------------------------------------
    nx = k - i1;
    b_y.set_size(nx);
    for (k = 0; k < nx; k++) {
      b_y[k] = clean_speech[i1 + k] * window[k];
    }

    fft(b_y, n_fft, b_x);
    nx = b_x.size(0);
    b_y.set_size(b_x.size(0));
    for (k = 0; k < nx; k++) {
      b_y[k] = rt_hypotd_snf(b_x[k].re, b_x[k].im);
    }

    i1 = b_y.size(0);
    clean_spec.set_size(b_y.size(0));
    for (k = 0; k < i1; k++) {
      clean_spec[k] = rt_powd_snf(b_y[k], 2.0);
    }

    nx = i2 - b_i;
    b_y.set_size(nx);
    for (i1 = 0; i1 < nx; i1++) {
      b_y[i1] = processed_speech[b_i + i1] * window[i1];
    }

    fft(b_y, n_fft, b_x);
    nx = b_x.size(0);
    b_y.set_size(b_x.size(0));
    for (k = 0; k < nx; k++) {
      b_y[k] = rt_hypotd_snf(b_x[k].re, b_x[k].im);
    }

    i1 = b_y.size(0);
    processed_spec.set_size(b_y.size(0));
    for (k = 0; k < i1; k++) {
      processed_spec[k] = rt_powd_snf(b_y[k], 2.0);
    }

    //  ----------------------------------------------------------
    //  (3) Compute Filterbank Output Energies (in dB scale)
    //  ----------------------------------------------------------
    for (b_i = 0; b_i < 25; b_i++) {
      b_y.set_size(loop_ub);
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_y[i1] = clean_spec[i1] * crit_filter[b_i + 25 * i1];
      }

      nx = b_y.size(0);
      if (b_y.size(0) == 0) {
        dBMax_processed = 0.0;
      } else {
        dBMax_processed = b_y[0];
        for (k = 2; k <= nx; k++) {
          dBMax_processed += b_y[k - 1];
        }
      }

      cent_freq[b_i] = dBMax_processed;
      b_y.set_size(b_loop_ub);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        b_y[i1] = processed_spec[i1] * crit_filter[b_i + 25 * i1];
      }

      nx = b_y.size(0);
      if (b_y.size(0) == 0) {
        dBMax_processed = 0.0;
      } else {
        dBMax_processed = b_y[0];
        for (k = 2; k <= nx; k++) {
          dBMax_processed += b_y[k - 1];
        }
      }

      absn = cent_freq[b_i];
      if (!(absn > 1.0E-10)) {
        absn = 1.0E-10;
      }

      cent_freq[b_i] = 10.0 * std::log10(absn);
      if (!(dBMax_processed > 1.0E-10)) {
        dBMax_processed = 1.0E-10;
      }

      bandwidth[b_i] = 10.0 * std::log10(dBMax_processed);
    }

    //  ----------------------------------------------------------
    //  (4) Compute Spectral Slope (dB[i+1]-dB[i])
    //  ----------------------------------------------------------
    for (i1 = 0; i1 < 24; i1++) {
      clean_slope[i1] = cent_freq[i1 + 1] - cent_freq[i1];
      processed_slope[i1] = bandwidth[i1 + 1] - bandwidth[i1];
    }

    //  ----------------------------------------------------------
    //  (5) Find the nearest peak locations in the spectra to
    //      each critical band.  If the slope is negative, we
    //      search to the left.  If positive, we search to the
    //      right.
    //  ----------------------------------------------------------
    //  ----------------------------------------------------------
    //   (6) Compute the WSS Measure for this frame.  This
    //       includes determination of the weighting function.
    //  ----------------------------------------------------------
    absn = maximum(cent_freq);
    dBMax_processed = maximum(bandwidth);

    //  The weights are calculated by averaging individual
    //  weighting factors from the clean and processed frame.
    //  These weights W_clean and W_processed should range
    //  from 0 to 1 and place more emphasis on spectral
    //  peaks and less emphasis on slope differences in spectral
    //  valleys.  This procedure is described on page 1280 of
    //  Klatt's 1982 ICASSP paper.
    for (b_i = 0; b_i < 24; b_i++) {
      //  find the peaks in the clean speech signal
      if (clean_slope[b_i] > 0.0) {
        //  search to the right
        nx = b_i;
        while ((nx + 1 < 25) && (clean_slope[nx] > 0.0)) {
          nx++;
        }

        clean_loc_peak[b_i] = cent_freq[nx - 1];
      } else {
        //  search to the left
        nx = b_i + 1;
        while ((nx > 0) && (clean_slope[nx - 1] <= 0.0)) {
          nx--;
        }

        clean_loc_peak[b_i] = cent_freq[nx];
      }

      //  find the peaks in the processed speech signal
      if (processed_slope[b_i] > 0.0) {
        //  search to the right
        nx = b_i;
        while ((nx + 1 < 25) && (processed_slope[nx] > 0.0)) {
          nx++;
        }

        x = bandwidth[nx - 1];
      } else {
        //  search to the left
        nx = b_i + 1;
        while ((nx > 0) && (processed_slope[nx - 1] <= 0.0)) {
          nx--;
        }

        x = bandwidth[nx];
      }

      clean_loc_peak[b_i] = (20.0 / ((absn + 20.0) - cent_freq[b_i]) * (1.0 /
        ((clean_loc_peak[b_i] + 1.0) - cent_freq[b_i])) + 20.0 /
        ((dBMax_processed + 20.0) - bandwidth[b_i]) * (1.0 / ((x + 1.0) -
        bandwidth[b_i]))) / 2.0;
    }

    for (k = 0; k < 24; k++) {
      x = clean_slope[k] - processed_slope[k];
      clean_slope[k] = x;
      processed_slope[k] = clean_loc_peak[k] * rt_powd_snf(x, 2.0);
    }

    dBMax_processed = processed_slope[0];

    //  this normalization is not part of Klatt's paper, but helps
    //  to normalize the measure.  Here we scale the measure by the
    //  sum of the weights.
    absn = clean_loc_peak[0];
    for (k = 0; k < 23; k++) {
      dBMax_processed += processed_slope[k + 1];
      absn += clean_loc_peak[k + 1];
    }

    distortion[frame_count] = dBMax_processed / absn;
    max_freq += skiprate;
  }
}

void composite_cpp(const coder::array<double, 1U> &cleanFile, const coder::array<
                   double, 1U> &enhancedFile, double Srate1, double Srate2, double
                   Csig_data[], int Csig_size[2], double Cbak_data[], int
                   Cbak_size[2], double Covl_data[], int Covl_size[2])
{
  int LLR_len;
  int len;
  int loop_ub;
  coder::array<double, 1U> data1;
  int i;
  coder::array<double, 1U> data2;
  coder::array<double, 2U> wss_dist_vec;
  double y;
  double wss_dist;
  coder::array<double, 2U> LLRs;
  double llr_mean;
  double winlength;
  double skiprate;
  double num_frames;
  double start;
  coder::array<double, 1U> window;
  double d;
  double pesq_mos_scores_data[2];
  int pesq_mos_scores_size[2];
  coder::array<double, 1U> clean_frame;
  coder::array<double, 1U> b_y;
  double y_data[2];
  if (!isInitialized_composite_cpp) {
    composite_cpp_initialize();
  }

  //  ----------------------------------------------------------------------
  //           Composite Objective Speech Quality Measure
  //
  //    This function implements the composite objective measure proposed in
  //    [1].
  //
  //    Usage:  [sig,bak,ovl]=composite(cleanFile.wav, enhancedFile.wav)
  //
  //          cleanFile.wav - clean input file in .wav format
  //          enhancedFile  - enhanced output file in .wav format
  //          sig           - predicted rating [1-5] of speech distortion
  //          bak           - predicted rating [1-5] of noise distortion
  //          ovl           - predicted rating [1-5] of overall quality
  //
  //        In addition to the above ratings (sig, bak, & ovl) it returns
  //        the individual values of the LLR, SNRseg, WSS and PESQ measures.
  //
  //   Example call:  [sig,bak,ovl] =composite('sp04.wav','enhanced.wav')
  //
  //
  //   References:
  //
  //      [1]   Hu, Y. and Loizou, P. (2006). Evaluation of objective measures
  //            for speech enhancement. Proc. Interspeech, Pittsburg, PA.
  //
  //    Authors: Yi Hu and Philipos C. Loizou
  //    (the LLR, SNRseg and WSS measures were based on Bryan Pellom and John
  //      Hansen's implementations)
  //
  //  Copyright (c) 2006 by Philipos C. Loizou
  //  $Revision: 0.0 $  $Date: 10/09/2006 $
  //  ----------------------------------------------------------------------
  //  if nargin~=2
  //      fprintf('USAGE: [sig,bak,ovl]=composite(cleanFile.wav, enhancedFile.wav)\n'); 
  //      fprintf('For more help, type: help composite\n\n');
  //      return;
  //  end
  //  [data1, Srate1]= audioread(cleanFile);
  //  [data2, Srate2]= audioread(enhancedFile);
  LLR_len = cleanFile.size(0);
  len = enhancedFile.size(0);
  if (LLR_len < len) {
    len = LLR_len;
  }

  if (1 > len) {
    loop_ub = 0;
  } else {
    loop_ub = len;
  }

  data1.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    data1[i] = cleanFile[i] + 2.2204460492503131E-16;
  }

  if (1 > len) {
    loop_ub = 0;
  } else {
    loop_ub = len;
  }

  data2.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    data2[i] = enhancedFile[i] + 2.2204460492503131E-16;
  }

  //  -- compute the WSS measure ---
  //
  wss(data1, data2, Srate1, wss_dist_vec);
  sort(wss_dist_vec);
  i = static_cast<int>(rt_roundd_snf(static_cast<double>(wss_dist_vec.size(1)) *
    0.95));
  if (1 > i) {
    i = 0;
  }

  if (i == 0) {
    y = 0.0;
  } else {
    y = wss_dist_vec[0];
    for (len = 2; len <= i; len++) {
      y += wss_dist_vec[len - 1];
    }
  }

  wss_dist = y / static_cast<double>(i);

  //  --- compute the LLR measure ---------
  //
  llr(data1, data2, Srate1, wss_dist_vec);
  LLRs.set_size(1, wss_dist_vec.size(1));
  loop_ub = wss_dist_vec.size(0) * wss_dist_vec.size(1);
  for (i = 0; i < loop_ub; i++) {
    LLRs[i] = wss_dist_vec[i];
  }

  sort(LLRs);
  LLR_len = static_cast<int>(rt_roundd_snf(static_cast<double>(wss_dist_vec.size
    (1)) * 0.95));
  if (1 > LLR_len) {
    i = 0;
  } else {
    i = LLR_len;
  }

  if (i == 0) {
    y = 0.0;
  } else {
    y = LLRs[0];
    for (len = 2; len <= i; len++) {
      y += LLRs[len - 1];
    }
  }

  llr_mean = y / static_cast<double>(i);

  //  --- compute the SNRseg ----------------
  //
  //  ----------------------------------------------------------------------
  //  ----------------------------------------------------------------------
  //  Check the length of the clean and processed speech.  Must be the same.
  //  ----------------------------------------------------------------------
  //  if (clean_length ~= processed_length)
  //    disp('Error: Both Speech Files must be same length.');
  //    return
  //  end
  //  ----------------------------------------------------------------------
  //  Scale both clean speech and processed speech to have same dynamic
  //  range.  Also remove DC component from each signal
  //  ----------------------------------------------------------------------
  // clean_speech     = clean_speech     - mean(clean_speech);
  // processed_speech = processed_speech - mean(processed_speech);
  // processed_speech = processed_speech.*(max(abs(clean_speech))/ max(abs(processed_speech))); 
  //  ----------------------------------------------------------------------
  //  Global Variables
  //  ----------------------------------------------------------------------
  winlength = rt_roundd_snf(30.0 * Srate1 / 1000.0);

  // 240;		   % window length in samples
  skiprate = std::floor(winlength / 4.0);

  //  window skip in samples
  //  minimum SNR in dB
  //  maximum SNR in dB
  //  ----------------------------------------------------------------------
  //  For each frame of input speech, calculate the Segmental SNR
  //  ----------------------------------------------------------------------
  num_frames = static_cast<double>(data1.size(0)) / skiprate - winlength /
    skiprate;

  //  number of frames
  start = 1.0;

  //  starting sample
  if (winlength < 1.0) {
    wss_dist_vec.set_size(1, 0);
  } else if (rtIsInf(winlength) && (1.0 == winlength)) {
    wss_dist_vec.set_size(1, 1);
    wss_dist_vec[0] = rtNaN;
  } else {
    loop_ub = static_cast<int>(winlength - 1.0);
    wss_dist_vec.set_size(1, (loop_ub + 1));
    for (i = 0; i <= loop_ub; i++) {
      wss_dist_vec[i] = static_cast<double>(i) + 1.0;
    }
  }

  window.set_size(wss_dist_vec.size(1));
  loop_ub = wss_dist_vec.size(1);
  for (i = 0; i < loop_ub; i++) {
    window[i] = 6.2831853071795862 * wss_dist_vec[i] / (winlength + 1.0);
  }

  LLR_len = window.size(0);
  for (len = 0; len < LLR_len; len++) {
    window[len] = std::cos(window[len]);
  }

  loop_ub = window.size(0);
  for (i = 0; i < loop_ub; i++) {
    window[i] = 0.5 * (1.0 - window[i]);
  }

  loop_ub = static_cast<int>(std::floor(num_frames));
  wss_dist_vec.set_size(1, loop_ub);
  for (i = 0; i < loop_ub; i++) {
    wss_dist_vec[i] = 0.0;
  }

  i = static_cast<int>(num_frames);
  for (int frame_count = 0; frame_count < i; frame_count++) {
    int i1;

    //  ----------------------------------------------------------
    //  (1) Get the Frames for the test and reference speech.
    //      Multiply by Hanning Window.
    //  ----------------------------------------------------------
    d = (start + winlength) - 1.0;
    if (start > d) {
      LLR_len = 0;
      len = 0;
      i1 = 1;
    } else {
      i1 = static_cast<int>(start);
      LLR_len = i1 - 1;
      len = static_cast<int>(d);
    }

    loop_ub = len - LLR_len;
    clean_frame.set_size(loop_ub);
    for (len = 0; len < loop_ub; len++) {
      clean_frame[len] = data1[LLR_len + len] * window[len];
    }

    //  ----------------------------------------------------------
    //  (2) Compute the Segmental SNR
    //  ----------------------------------------------------------
    b_y.set_size(clean_frame.size(0));
    LLR_len = clean_frame.size(0);
    for (len = 0; len < LLR_len; len++) {
      b_y[len] = rt_powd_snf(clean_frame[len], 2.0);
    }

    LLR_len = b_y.size(0);
    if (b_y.size(0) == 0) {
      y = 0.0;
    } else {
      y = b_y[0];
      for (len = 2; len <= LLR_len; len++) {
        y += b_y[len - 1];
      }
    }

    loop_ub = clean_frame.size(0);
    for (LLR_len = 0; LLR_len < loop_ub; LLR_len++) {
      clean_frame[LLR_len] = clean_frame[LLR_len] - data2[(i1 + LLR_len) - 1] *
        window[LLR_len];
    }

    b_y.set_size(clean_frame.size(0));
    LLR_len = clean_frame.size(0);
    for (len = 0; len < LLR_len; len++) {
      b_y[len] = rt_powd_snf(clean_frame[len], 2.0);
    }

    LLR_len = b_y.size(0);
    if (b_y.size(0) == 0) {
      num_frames = 0.0;
    } else {
      num_frames = b_y[0];
      for (len = 2; len <= LLR_len; len++) {
        num_frames += b_y[len - 1];
      }
    }

    d = 10.0 * std::log10(y / (num_frames + 2.2204460492503131E-16) +
                          2.2204460492503131E-16);
    wss_dist_vec[frame_count] = d;
    if (!(d > -10.0)) {
      d = -10.0;
      wss_dist_vec[frame_count] = -10.0;
    }

    if (!(d < 35.0)) {
      wss_dist_vec[frame_count] = 35.0;
    }

    start += skiprate;
  }

  LLR_len = wss_dist_vec.size(1);
  if (wss_dist_vec.size(1) == 0) {
    y = 0.0;
  } else {
    y = wss_dist_vec[0];
    for (len = 2; len <= LLR_len; len++) {
      y += wss_dist_vec[len - 1];
    }
  }

  //  -- compute the pesq ----
  //
  //  if     Srate1==8000,    mode='nb';
  //  elseif Srate1 == 16000, mode='wb';
  //  else,
  //       error ('Sampling freq in PESQ needs to be 8 kHz or 16 kHz');
  //  end
  pesq_original_cpp(cleanFile, Srate1, enhancedFile, pesq_mos_scores_data,
                    pesq_mos_scores_size);
  if (pesq_mos_scores_size[1] == 2) {
    d = pesq_mos_scores_data[0];
    pesq_mos_scores_size[0] = 1;
    pesq_mos_scores_size[1] = 1;
    pesq_mos_scores_data[0] = d;

    //  take the raw PESQ value instead of the
    //  MOS-mapped value (this composite
    //  measure was only validated with the raw
    //  PESQ value)
  }

  //  --- now compute the composite measures ------------------
  //
  d = 3.093 - 1.029 * llr_mean;
  num_frames = 0.009 * wss_dist;
  LLR_len = pesq_mos_scores_size[0] * pesq_mos_scores_size[1];
  for (i = 0; i < LLR_len; i++) {
    y_data[i] = (d + 0.603 * pesq_mos_scores_data[i]) - num_frames;
  }

  Csig_size[1] = static_cast<signed char>(pesq_mos_scores_size[1]);
  if (0 <= static_cast<signed char>(pesq_mos_scores_size[1]) - 1) {
    if ((1.0 > y_data[0]) || rtIsNaN(y_data[0])) {
      Csig_data[0] = 1.0;
    } else {
      Csig_data[0] = y_data[0];
    }
  }

  loop_ub = Csig_size[1];
  if (0 <= loop_ub - 1) {
    std::memcpy(&y_data[0], &Csig_data[0], loop_ub * sizeof(double));
  }

  Csig_size[0] = 1;
  if (0 <= Csig_size[1] - 1) {
    if ((5.0 < y_data[0]) || rtIsNaN(y_data[0])) {
      Csig_data[0] = 5.0;
    } else {
      Csig_data[0] = y_data[0];
    }
  }

  //  limit values to [1, 5]
  d = 0.007 * wss_dist;
  num_frames = 0.063 * (y / static_cast<double>(wss_dist_vec.size(1)));
  for (i = 0; i < LLR_len; i++) {
    y_data[i] = ((0.478 * pesq_mos_scores_data[i] + 1.634) - d) + num_frames;
  }

  Cbak_size[1] = static_cast<signed char>(pesq_mos_scores_size[1]);
  if (0 <= static_cast<signed char>(pesq_mos_scores_size[1]) - 1) {
    if ((1.0 > y_data[0]) || rtIsNaN(y_data[0])) {
      Cbak_data[0] = 1.0;
    } else {
      Cbak_data[0] = y_data[0];
    }
  }

  loop_ub = Cbak_size[1];
  if (0 <= loop_ub - 1) {
    std::memcpy(&y_data[0], &Cbak_data[0], loop_ub * sizeof(double));
  }

  Cbak_size[0] = 1;
  if (0 <= Cbak_size[1] - 1) {
    if ((5.0 < y_data[0]) || rtIsNaN(y_data[0])) {
      Cbak_data[0] = 5.0;
    } else {
      Cbak_data[0] = y_data[0];
    }
  }

  //  limit values to [1, 5]
  d = 0.512 * llr_mean;
  num_frames = 0.007 * wss_dist;
  loop_ub = LLR_len - 1;
  for (i = 0; i <= loop_ub; i++) {
    pesq_mos_scores_data[i] = ((0.805 * pesq_mos_scores_data[i] + 1.594) - d) -
      num_frames;
  }

  i = static_cast<signed char>(pesq_mos_scores_size[1]) - 1;
  if (0 <= i) {
    if ((1.0 > pesq_mos_scores_data[0]) || rtIsNaN(pesq_mos_scores_data[0])) {
      Covl_data[0] = 1.0;
    } else {
      Covl_data[0] = pesq_mos_scores_data[0];
    }
  }

  loop_ub = static_cast<signed char>(pesq_mos_scores_size[1]);
  if (0 <= loop_ub - 1) {
    std::memcpy(&y_data[0], &Covl_data[0], loop_ub * sizeof(double));
  }

  Covl_size[0] = 1;
  Covl_size[1] = static_cast<signed char>(pesq_mos_scores_size[1]);
  if (0 <= i) {
    if ((5.0 < y_data[0]) || rtIsNaN(y_data[0])) {
      Covl_data[0] = 5.0;
    } else {
      Covl_data[0] = y_data[0];
    }
  }

  //  limit values to [1, 5]
  //  fprintf('\n LLR=%f   SNRseg=%f   WSS=%f   PESQ=%f\n',llr_mean,segSNR,wss_dist,pesq_mos); 
}

// End of code generation (composite_cpp.cpp)
