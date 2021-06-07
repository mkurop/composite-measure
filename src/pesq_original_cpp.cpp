//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  pesq_original_cpp.cpp
//
//  Code generation for function 'pesq_original_cpp'
//


// Include files
#include <pesq_original_cpp.h>
#include <FFTImplementationCallback.h>
#include <abs.h>
#include <composite_cpp.h>
#include <composite_cpp_rtwutil.h>
#include <fft.h>
#include <filter.h>
#include <find.h>
#include <gencoswin.h>
#include <hann.h>
#include <ifft.h>
#include <interp1.h>
#include <log2.h>
#include <minOrMax.h>
#include <relop.h>
#include <rt_defines.h>
#include <rt_nonfinite.h>
#include <sprintf.h>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <math.h>

// Type Definitions
struct emxArray_real_T_12x5
{
  double data[60];
  int size[2];
};

struct emxArray_real_T_1x49
{
  double data[49];
  int size[2];
};

struct emxArray_real_T_1x1024
{
  double data[1024];
  int size[2];
};

struct c_struct_T
{
  double CALIBRATE;
  double Nfmax;
  double MAXNUTTERANCES;
  double MINUTTLENGTH;
  double WHOLE_SIGNAL;
  double UttSearch_Start[50];
  double UttSearch_End[50];
  double Utt_DelayEst[50];
  double Utt_Delay[50];
  double Utt_DelayConf[50];
  double Utt_Start[50];
  double Utt_End[50];
  double DATAPADDING_MSEC;
  double SEARCHBUFFER;
  double MINSPEECHLGTH;
  double JOINSPEECHLGTH;
  double Crude_DelayEst;
  double Crude_DelayConf;
  double Nutterances;
  double Largest_uttsize;
  double Best_DC1;
  double Best_DC2;
  double Best_ED1;
  double Best_D1;
  double Best_ED2;
  double Best_D2;
  double Best_BP;
  double WB_InIIR_Nsos;
  double WB_InIIR_Hsos[5];
  double Downsample;
  double Nb;
  double Sl;
  double Sp;
  double Fs;
  emxArray_real_T_12x5 InIIR_Hsos;
  double InIIR_Nsos;
  double Align_Nfft;
  emxArray_real_T_1x49 nr_of_hz_bands_per_bark_band;
  emxArray_real_T_1x49 centre_of_band_bark;
  emxArray_real_T_1x49 centre_of_band_hz;
  emxArray_real_T_1x49 width_of_band_bark;
  emxArray_real_T_1x49 width_of_band_hz;
  emxArray_real_T_1x49 pow_dens_correction_factor;
  emxArray_real_T_1x49 abs_thresh_power;
  emxArray_real_T_1x1024 Window;
  double TARGET_AVG_POWER;
  double Plot_Frame;
  double c_NUMBER_OF_PSQM_FRAMES_PER_SYL;
};

struct struct_T
{
  double CALIBRATE;
  double Nfmax;
  double MAXNUTTERANCES;
  double MINUTTLENGTH;
  double WHOLE_SIGNAL;
  double UttSearch_Start[50];
  double UttSearch_End[50];
  double Utt_DelayEst[50];
  double Utt_Delay[50];
  double Utt_DelayConf[50];
  double Utt_Start[50];
  double Utt_End[50];
  double DATAPADDING_MSEC;
  double SEARCHBUFFER;
  double MINSPEECHLGTH;
  double JOINSPEECHLGTH;
  double Crude_DelayEst;
  double Crude_DelayConf;
  double Nutterances;
  double Largest_uttsize;
  double Best_DC1;
  double Best_DC2;
  double Best_ED1;
  double Best_D1;
  double Best_ED2;
  double Best_D2;
  double Best_BP;
  double WB_InIIR_Nsos;
  double WB_InIIR_Hsos[5];
  double Downsample;
  double Nb;
  double Sl;
  double Sp;
  double Fs;
  emxArray_real_T_12x5 InIIR_Hsos;
  double InIIR_Nsos;
  double Align_Nfft;
  emxArray_real_T_1x49 nr_of_hz_bands_per_bark_band;
  emxArray_real_T_1x49 centre_of_band_bark;
  emxArray_real_T_1x49 centre_of_band_hz;
  emxArray_real_T_1x49 width_of_band_bark;
  emxArray_real_T_1x49 width_of_band_hz;
  emxArray_real_T_1x49 pow_dens_correction_factor;
  emxArray_real_T_1x49 abs_thresh_power;
  emxArray_real_T_1x1024 Window;
  double TARGET_AVG_POWER;
};

struct b_struct_T
{
  double CALIBRATE;
  double Nfmax;
  double MAXNUTTERANCES;
  double MINUTTLENGTH;
  double WHOLE_SIGNAL;
  double UttSearch_Start[50];
  double UttSearch_End[50];
  double Utt_DelayEst[50];
  double Utt_Delay[50];
  double Utt_DelayConf[50];
  double Utt_Start[50];
  double Utt_End[50];
  double DATAPADDING_MSEC;
  double SEARCHBUFFER;
  double MINSPEECHLGTH;
  double JOINSPEECHLGTH;
  double Crude_DelayEst;
  double Crude_DelayConf;
  double Nutterances;
  double Largest_uttsize;
  double Best_DC1;
  double Best_DC2;
  double Best_ED1;
  double Best_D1;
  double Best_ED2;
  double Best_D2;
  double Best_BP;
  double WB_InIIR_Nsos;
  double WB_InIIR_Hsos[5];
  double Downsample;
  double Nb;
  double Sl;
  double Sp;
  double Fs;
  emxArray_real_T_12x5 InIIR_Hsos;
  double InIIR_Nsos;
  double Align_Nfft;
  emxArray_real_T_1x49 nr_of_hz_bands_per_bark_band;
  emxArray_real_T_1x49 centre_of_band_bark;
  emxArray_real_T_1x49 centre_of_band_hz;
  emxArray_real_T_1x49 width_of_band_bark;
  emxArray_real_T_1x49 width_of_band_hz;
  emxArray_real_T_1x49 pow_dens_correction_factor;
  emxArray_real_T_1x49 abs_thresh_power;
  emxArray_real_T_1x1024 Window;
};

// Function Declarations
static void DC_block(const coder::array<double, 2U> &data, double Nsamples,
                     const struct_T *glb, coder::array<double, 2U> &mod_data);
static double Lpq_weight(double start_frame, double stop_frame, const coder::
  array<double, 2U> &frame_disturbance, const coder::array<double, 2U>
  &time_weight);
static void apply_VAD(const coder::array<double, 2U> &data, double Nsamples,
                      const struct_T *glb, coder::array<double, 2U> &VAD, coder::
                      array<double, 2U> &logVAD);
static void apply_filter(const coder::array<double, 2U> &data, double
  data_Nsamples, const struct_T *glb, coder::array<double, 2U> &align_filtered);
static void apply_filters(const coder::array<double, 2U> &data, const struct_T
  *glb, coder::array<double, 2U> &mod_data);
static void apply_filters_WB(const coder::array<double, 2U> &data, const
  struct_T *glb, coder::array<double, 2U> &mod_data);
static void b_apply_filter(const coder::array<double, 2U> &data, double
  data_Nsamples, const struct_T *glb, coder::array<double, 2U> &align_filtered);
static void b_fix_power_level(const coder::array<double, 2U> &data, double
  data_Nsamples, double maxNsamples, struct_T *glb, coder::array<double, 2U>
  &mod_data);
static void compute_delay(double stop_sample, double search_range, const coder::
  array<double, 2U> &time_series1, const coder::array<double, 2U> &time_series2,
  double *best_delay, double *max_correlation);
static void crude_align(const coder::array<double, 2U> &ref_logVAD, double
  ref_Nsamples, const coder::array<double, 2U> &deg_logVAD, double deg_Nsamples,
  double Utt_id, struct_T *glb);
static void fix_power_level(const coder::array<double, 2U> &data, double
  data_Nsamples, double maxNsamples, const b_struct_T *glb, coder::array<double,
  2U> &mod_data, struct_T *b_glb);
static void freq_resp_compensation(double number_of_frames, const coder::array<
  double, 2U> &pitch_pow_dens_ref, const double avg_pitch_pow_dens_ref_data[],
  const double avg_pitch_pow_dens_deg_data[], const c_struct_T *glb, coder::
  array<double, 2U> &mod_pitch_pow_dens_ref);
static void freq_warping(const double hz_spectrum_data[], const c_struct_T *glb,
  double pitch_pow_dens_data[], int pitch_pow_dens_size[2]);
static void input_filter(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &deg_data, double deg_Nsamples,
  const struct_T *glb, coder::array<double, 2U> &mod_ref_data, coder::array<
  double, 2U> &mod_deg_data);
static void intensity_warping_of(double frame, const coder::array<double, 2U>
  &pitch_pow_dens, const c_struct_T *glb, double loudness_dens_data[], int
  loudness_dens_size[2]);
static void multiply_with_asymmetry_factor(const double disturbance_dens_data[],
  const int disturbance_dens_size[2], double frame, const coder::array<double,
  2U> &pitch_pow_dens_ref, const coder::array<double, 2U> &pitch_pow_dens_deg,
  const c_struct_T *glb, double mod_disturbance_dens_data[], int
  mod_disturbance_dens_size[2]);
static void pesq_psychoacoustic_model(const coder::array<double, 2U> &ref_data,
  double ref_Nsamples, const coder::array<double, 2U> &deg_data, double
  deg_Nsamples, const struct_T *glb, double *pesq_mos, c_struct_T *b_glb);
static double pseudo_Lp(const double x_data[], double p, const c_struct_T *glb);
static double rt_atan2d_snf(double u0, double u1);
static double rt_remd_snf(double u0, double u1);
static void short_term_fft(double Nf, const coder::array<double, 2U> &data,
  const double Whanning_data[], double start_sample, double hz_spectrum_data[],
  int hz_spectrum_size[2]);
static void split_align(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &ref_logVAD, const coder::array<
  double, 2U> &deg_data, double deg_Nsamples, const coder::array<double, 2U>
  &deg_logVAD, double Utt_Start_l, double Utt_SpeechStart, double Utt_SpeechEnd,
  double Utt_End_l, double Utt_DelayEst_l, double Utt_DelayConf_l, struct_T *glb);
static void time_align(const coder::array<double, 2U> &ref_data, const coder::
  array<double, 2U> &deg_data, double deg_Nsamples, double Utt_id, struct_T *glb);
static void time_avg_audible_of(double number_of_frames, const coder::array<
  double, 2U> &silent, const coder::array<double, 2U> &pitch_pow_dens, double
  total_number_of_frames, const c_struct_T *glb, double avg_pitch_pow_dens_data[],
  int avg_pitch_pow_dens_size[2]);
static double total_audible(double frame, const coder::array<double, 2U>
  &pitch_pow_dens, double factor, const c_struct_T *glb);
static void utterance_locate(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &ref_VAD, const coder::array<
  double, 2U> &ref_logVAD, const coder::array<double, 2U> &deg_data, double
  deg_Nsamples, const coder::array<double, 2U> &deg_logVAD, struct_T *glb);

// Function Definitions
static void DC_block(const coder::array<double, 2U> &data, double Nsamples,
                     const struct_T *glb, coder::array<double, 2U> &mod_data)
{
  double ofs;
  int loop_ub;
  int i;
  double d;
  int vlen_tmp;
  double facc;
  int k;
  int i1;
  int tmp_data[64];
  coder::array<double, 2U> y;
  double mod_data_data[64];
  int i2;
  coder::array<double, 2U> b_mod_data;

  //      global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER
  ofs = 75.0 * glb->Downsample;
  mod_data.set_size(1, data.size(1));
  loop_ub = data.size(0) * data.size(1);
  for (i = 0; i < loop_ub; i++) {
    mod_data[i] = data[i];
  }

  // compute dc component, it is a little weird
  d = Nsamples - ofs;
  if (ofs + 1.0 > d) {
    i = -1;
    vlen_tmp = -1;
  } else {
    i = static_cast<int>(ofs + 1.0) - 2;
    vlen_tmp = static_cast<int>(d) - 1;
  }

  vlen_tmp -= i;
  if (vlen_tmp == 0) {
    facc = 0.0;
  } else {
    facc = data[i + 1];
    for (k = 2; k <= vlen_tmp; k++) {
      facc += data[i + k];
    }
  }

  facc /= Nsamples;
  if (ofs + 1.0 > d) {
    vlen_tmp = 0;
    i1 = 0;
    i = 1;
  } else {
    i = static_cast<int>(ofs + 1.0);
    vlen_tmp = i - 1;
    i1 = static_cast<int>(d);
  }

  loop_ub = i1 - vlen_tmp;
  for (i1 = 0; i1 < loop_ub; i1++) {
    mod_data[(i + i1) - 1] = data[vlen_tmp + i1] - facc;
  }

  loop_ub = static_cast<int>(std::floor(glb->Downsample - 1.0));
  vlen_tmp = loop_ub + 1;
  for (i = 0; i <= loop_ub; i++) {
    tmp_data[i] = static_cast<int>(ofs + (static_cast<double>(i) + 1.0));
  }

  if (rtIsNaN(glb->Downsample - 1.0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else if (glb->Downsample - 1.0 < 0.0) {
    y.set_size(1, 0);
  } else if (rtIsInf(glb->Downsample - 1.0) && (0.0 == glb->Downsample - 1.0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else {
    facc = glb->Downsample - 1.0;
    y.set_size(1, (static_cast<int>(std::floor(facc)) + 1));
    loop_ub = static_cast<int>(std::floor(facc));
    for (i = 0; i <= loop_ub; i++) {
      y[i] = i;
    }
  }

  for (i = 0; i < vlen_tmp; i++) {
    mod_data_data[i] = mod_data[tmp_data[i] - 1] * (y[i] + 0.5) /
      glb->Downsample;
  }

  for (i = 0; i < vlen_tmp; i++) {
    mod_data[tmp_data[i] - 1] = mod_data_data[i];
  }

  facc = (d - glb->Downsample) + 1.0;
  if (facc > d) {
    i = 0;
    vlen_tmp = 1;
    i1 = -1;
    k = 1;
    i2 = 1;
  } else {
    i = static_cast<int>(d) - 1;
    vlen_tmp = -1;
    i1 = static_cast<int>(facc) - 1;
    k = static_cast<int>(d);
    i2 = -1;
  }

  loop_ub = div_s32_floor(i1 - i, vlen_tmp);
  b_mod_data.set_size(1, (loop_ub + 1));
  for (i1 = 0; i1 <= loop_ub; i1++) {
    b_mod_data[i1] = mod_data[i + vlen_tmp * i1] * (y[i1] + 0.5) /
      glb->Downsample;
  }

  loop_ub = b_mod_data.size(1);
  for (i = 0; i < loop_ub; i++) {
    mod_data[(k + i2 * i) - 1] = b_mod_data[i];
  }
}

static double Lpq_weight(double start_frame, double stop_frame, const coder::
  array<double, 2U> &frame_disturbance, const coder::array<double, 2U>
  &time_weight)
{
  double result_time;
  double total_time_weight_time;
  int i;
  double result_syllable;

  //  fid= fopen( 'tmp_mat.txt', 'wt');
  //  fprintf( fid, 'time_weight\n');
  //  fprintf( fid, '%f\n', time_weight);
  //  fprintf( fid, 'frame_disturbance:\n');
  //  fprintf( fid, '%f\n', frame_disturbance);
  //  fprintf( fid, 'frame_disturbance_asym_add\n');
  //  fprintf( fid, '%f\n', frame_disturbance_asym_add);
  //  fclose( fid);
  //      global glb.NUMBER_OF_PSQM_FRAMES_PER_SYLLABE
  //  fid= fopen( 'tmp_mat1.txt', 'at');
  //  fprintf( 'result_time:\n');
  result_time = 0.0;
  total_time_weight_time = 0.0;

  //  fprintf( 'start/end frame: %d/%d\n', start_frame, stop_frame);
  i = static_cast<int>((stop_frame + (10.0 - start_frame)) / 10.0);
  for (int start_frame_of_syllable = 0; start_frame_of_syllable < i;
       start_frame_of_syllable++) {
    double b_start_frame_of_syllable;
    double count_syllable;
    int i1;
    double b_frame;
    b_start_frame_of_syllable = start_frame + static_cast<double>
      (start_frame_of_syllable) * 10.0;
    result_syllable = 0.0;
    count_syllable = 0.0;
    i1 = static_cast<int>(((b_start_frame_of_syllable + 20.0) - 1.0) + (1.0 -
      b_start_frame_of_syllable));
    for (int frame = 0; frame < i1; frame++) {
      b_frame = b_start_frame_of_syllable + static_cast<double>(frame);
      if (b_frame <= stop_frame) {
        //              if (start_frame_of_syllable== 101)
        //                  fprintf( fid, '%f\n', h);
        //              end
        result_syllable += rt_powd_snf(frame_disturbance[static_cast<int>
          (b_frame + 1.0) - 1], 6.0);
      }

      count_syllable++;
    }

    b_frame = time_weight[static_cast<int>((b_start_frame_of_syllable + 1.0) -
      start_frame) - 1];
    b_start_frame_of_syllable = b_frame * rt_powd_snf(result_syllable /
      count_syllable, 0.16666666666666666);
    result_time += b_start_frame_of_syllable * b_start_frame_of_syllable;
    total_time_weight_time += b_frame * b_frame;

    //      fprintf( fid, '%f\n', result_time);
  }

  //  fclose (fid);
  //  fprintf( 'total_time_weight_time is %f\n', total_time_weight_time);
  result_time /= total_time_weight_time;
  return std::sqrt(result_time);
}

static void apply_VAD(const coder::array<double, 2U> &data, double Nsamples,
                      const struct_T *glb, coder::array<double, 2U> &VAD, coder::
                      array<double, 2U> &logVAD)
{
  double Nwindows;
  int i;
  int count;
  double d;
  int vlen;
  double StDNoise;
  int i1;
  int k;
  int iteration;
  double LevelThresh;
  double LevelMin;
  coder::array<double, 2U> y;
  int end;
  int start;
  coder::array<double, 2U> VAD_lessthan_LevelThresh;
  coder::array<double, 2U> VAD_greaterthan_LevelThresh;
  double LevelNoise;
  double LevelSig;
  coder::array<double, 2U> a;
  int loop_ub;
  coder::array<int, 2U> r;
  coder::array<int, 2U> r1;
  coder::array<double, 1U> b_VAD;
  coder::array<boolean_T, 2U> c_VAD;

  //  -------------
  //      global glb.Downsample glb.MINSPEECHLGTH glb.JOINSPEECHLGTH
  Nwindows = std::floor(Nsamples / glb->Downsample);

  // number of 4ms window
  i = static_cast<int>(Nwindows);
  VAD.set_size(1, i);
  for (count = 0; count < i; count++) {
    d = ((static_cast<double>(count) + 1.0) - 1.0) * glb->Downsample + 1.0;
    StDNoise = (static_cast<double>(count) + 1.0) * glb->Downsample;
    if (d > StDNoise) {
      i1 = -1;
      iteration = 0;
    } else {
      i1 = static_cast<int>(d) - 2;
      iteration = static_cast<int>(StDNoise);
    }

    iteration = (iteration - i1) - 1;
    y.set_size(1, iteration);
    for (k = 0; k < iteration; k++) {
      y[k] = rt_powd_snf(data[(i1 + k) + 1], 2.0);
    }

    vlen = y.size(1);
    if (y.size(1) == 0) {
      StDNoise = 0.0;
    } else {
      StDNoise = y[0];
      for (k = 2; k <= vlen; k++) {
        StDNoise += y[k - 1];
      }
    }

    VAD[count] = StDNoise / glb->Downsample;
  }

  // VAD is the power of each 4ms window
  vlen = VAD.size(1);
  if (VAD.size(1) == 0) {
    StDNoise = 0.0;
  } else {
    StDNoise = VAD[0];
    for (k = 2; k <= vlen; k++) {
      StDNoise += VAD[k - 1];
    }
  }

  LevelThresh = StDNoise / Nwindows;

  // LevelThresh is set to mean value of VAD
  LevelMin = b_maximum(VAD);
  if (LevelMin > 0.0) {
    LevelMin *= 0.0001;
  } else {
    LevelMin = 1.0;
  }

  // fprintf( 1, 'LevelMin is %f\n', LevelMin);
  end = VAD.size(1);
  for (start = 0; start < end; start++) {
    if (VAD[start] < LevelMin) {
      VAD[start] = LevelMin;
    }
  }

  end = VAD.size(1) - 1;
  for (iteration = 0; iteration < 12; iteration++) {
    StDNoise = 0.0;
    vlen = 0;
    for (start = 0; start <= end; start++) {
      if (VAD[start] <= LevelThresh) {
        vlen++;
      }
    }

    VAD_lessthan_LevelThresh.set_size(1, vlen);
    vlen = 0;
    for (start = 0; start <= end; start++) {
      d = VAD[start];
      if (d <= LevelThresh) {
        VAD_lessthan_LevelThresh[vlen] = d;
        vlen++;
      }
    }

    i1 = VAD_lessthan_LevelThresh.size(1);
    if (VAD_lessthan_LevelThresh.size(1) == 0) {
      LevelNoise = 0.0;
    } else {
      LevelNoise = VAD_lessthan_LevelThresh[0];
      for (k = 2; k <= i1; k++) {
        LevelNoise += VAD_lessthan_LevelThresh[k - 1];
      }
    }

    if (VAD_lessthan_LevelThresh.size(1) > 0) {
      LevelNoise /= static_cast<double>(VAD_lessthan_LevelThresh.size(1));
      a.set_size(1, VAD_lessthan_LevelThresh.size(1));
      loop_ub = VAD_lessthan_LevelThresh.size(0) * VAD_lessthan_LevelThresh.size
        (1);
      for (i1 = 0; i1 < loop_ub; i1++) {
        a[i1] = VAD_lessthan_LevelThresh[i1] - LevelNoise;
      }

      y.set_size(1, a.size(1));
      vlen = a.size(1);
      for (k = 0; k < vlen; k++) {
        y[k] = rt_powd_snf(a[k], 2.0);
      }

      vlen = y.size(1);
      StDNoise = y[0];
      for (k = 2; k <= vlen; k++) {
        StDNoise += y[k - 1];
      }

      StDNoise = std::sqrt(StDNoise / static_cast<double>
                           (VAD_lessthan_LevelThresh.size(1)));
    }

    LevelThresh = 1.001 * (LevelNoise + 2.0 * StDNoise);
  }

  // fprintf( 1, 'LevelThresh is %f\n', LevelThresh);
  end = VAD.size(1) - 1;
  vlen = 0;
  for (start = 0; start <= end; start++) {
    if (VAD[start] > LevelThresh) {
      vlen++;
    }
  }

  VAD_greaterthan_LevelThresh.set_size(1, vlen);
  vlen = 0;
  for (start = 0; start <= end; start++) {
    d = VAD[start];
    if (d > LevelThresh) {
      VAD_greaterthan_LevelThresh[vlen] = d;
      vlen++;
    }
  }

  vlen = VAD_greaterthan_LevelThresh.size(1);
  if (VAD_greaterthan_LevelThresh.size(1) == 0) {
    LevelSig = 0.0;
  } else {
    LevelSig = VAD_greaterthan_LevelThresh[0];
    for (k = 2; k <= vlen; k++) {
      LevelSig += VAD_greaterthan_LevelThresh[k - 1];
    }
  }

  end = VAD.size(1) - 1;
  vlen = 0;
  for (start = 0; start <= end; start++) {
    if (VAD[start] <= LevelThresh) {
      vlen++;
    }
  }

  r.set_size(1, vlen);
  vlen = 0;
  for (start = 0; start <= end; start++) {
    if (VAD[start] <= LevelThresh) {
      r[vlen] = start + 1;
      vlen++;
    }
  }

  if (VAD_greaterthan_LevelThresh.size(1) > 0) {
    LevelSig /= static_cast<double>(VAD_greaterthan_LevelThresh.size(1));
  } else {
    LevelThresh = -1.0;
  }

  // fprintf( 1, 'LevelSig is %f\n', LevelSig);
  if (VAD_greaterthan_LevelThresh.size(1) < Nwindows) {
    y.set_size(1, r.size(1));
    loop_ub = r.size(0) * r.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      y[i1] = VAD[r[i1] - 1];
    }

    vlen = r.size(1);
    if (r.size(1) == 0) {
      StDNoise = 0.0;
    } else {
      StDNoise = VAD[r[0] - 1];
      for (k = 2; k <= vlen; k++) {
        StDNoise += y[k - 1];
      }
    }

    LevelNoise = StDNoise / (Nwindows - static_cast<double>
      (VAD_greaterthan_LevelThresh.size(1)));
  } else {
    LevelNoise = 1.0;
  }

  // fprintf( 1, 'LevelNoise is %f\n', LevelNoise);
  end = VAD.size(1) - 1;
  vlen = 0;
  for (start = 0; start <= end; start++) {
    if (VAD[start] <= LevelThresh) {
      vlen++;
    }
  }

  r1.set_size(1, vlen);
  vlen = 0;
  for (start = 0; start <= end; start++) {
    if (VAD[start] <= LevelThresh) {
      r1[vlen] = start + 1;
      vlen++;
    }
  }

  loop_ub = r1.size(0) * r1.size(1);
  b_VAD.set_size(loop_ub);
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_VAD[i1] = -VAD[r1[i1] - 1];
  }

  loop_ub = b_VAD.size(0);
  for (i1 = 0; i1 < loop_ub; i1++) {
    VAD[r1[i1] - 1] = b_VAD[i1];
  }

  VAD[0] = -LevelMin;
  end = i - 1;
  VAD[end] = -LevelMin;
  start = -1;
  i = static_cast<int>(Nwindows + -1.0);
  for (count = 0; count < i; count++) {
    d = VAD[count + 1];
    if ((d > 0.0) && (VAD[count] <= 0.0)) {
      start = count + 1;
    }

    if ((d <= 0.0) && (VAD[count] > 0.0) && ((count - start) + 1 <= 4)) {
      if (start + 1 > count + 1) {
        i1 = 0;
        iteration = -1;
        vlen = 0;
      } else {
        i1 = start;
        iteration = count;
        vlen = start;
      }

      loop_ub = iteration - i1;
      y.set_size(1, (loop_ub + 1));
      for (iteration = 0; iteration <= loop_ub; iteration++) {
        y[iteration] = -VAD[i1 + iteration];
      }

      loop_ub = y.size(1);
      for (i1 = 0; i1 < loop_ub; i1++) {
        VAD[vlen + i1] = y[i1];
      }
    }
  }

  // to make sure finish- start is more than 4
  if (LevelSig >= LevelNoise * 1000.0) {
    for (count = 0; count < i; count++) {
      d = VAD[count + 1];
      if ((d > 0.0) && (VAD[count] <= 0.0)) {
        start = count + 1;
      }

      if ((d <= 0.0) && (VAD[count] > 0.0)) {
        if (start + 1 > count + 1) {
          i1 = -1;
          iteration = -1;
        } else {
          i1 = start - 1;
          iteration = count;
        }

        loop_ub = iteration - i1;
        y.set_size(1, loop_ub);
        for (iteration = 0; iteration < loop_ub; iteration++) {
          y[iteration] = VAD[(i1 + iteration) + 1];
        }

        if (loop_ub == 0) {
          StDNoise = 0.0;
        } else {
          StDNoise = VAD[i1 + 1];
          for (k = 2; k <= loop_ub; k++) {
            StDNoise += y[k - 1];
          }
        }

        if (StDNoise < 3.0 * LevelThresh * static_cast<double>((count - start) +
             1)) {
          if (start + 1 > count + 1) {
            i1 = 0;
            iteration = -1;
            vlen = 0;
          } else {
            i1 = start;
            iteration = count;
            vlen = start;
          }

          loop_ub = iteration - i1;
          y.set_size(1, (loop_ub + 1));
          for (iteration = 0; iteration <= loop_ub; iteration++) {
            y[iteration] = -VAD[i1 + iteration];
          }

          loop_ub = y.size(1);
          for (i1 = 0; i1 < loop_ub; i1++) {
            VAD[vlen + i1] = y[i1];
          }
        }
      }
    }
  }

  vlen = -2;
  for (count = 0; count < i; count++) {
    if ((VAD[count + 1] > 0.0) && (VAD[count] <= 0.0) && (vlen + 2 > 0) &&
        (count - vlen <= 50)) {
      if (vlen + 2 > count + 1) {
        i1 = -1;
        iteration = -1;
      } else {
        i1 = vlen;
        iteration = count;
      }

      loop_ub = iteration - i1;
      for (iteration = 0; iteration < loop_ub; iteration++) {
        VAD[(i1 + iteration) + 1] = LevelMin;
      }
    }

    if ((VAD[count + 1] <= 0.0) && (VAD[count] > 0.0)) {
      vlen = count;
    }
  }

  start = 0;
  for (count = 0; count < i; count++) {
    if ((VAD[count + 1] > 0.0) && (VAD[count] <= 0.0)) {
      start = count + 2;
    }
  }

  if (start == 0) {
    vlen = VAD.size(1);
    y.set_size(1, VAD.size(1));
    for (k = 0; k < vlen; k++) {
      y[k] = std::abs(VAD[k]);
    }

    VAD.set_size(1, y.size(1));
    loop_ub = y.size(0) * y.size(1);
    for (i = 0; i < loop_ub; i++) {
      VAD[i] = y[i];
    }

    VAD[0] = -LevelMin;
    VAD[end] = -LevelMin;
  }

  count = 3;
  while (count + 1 < Nwindows - 1.0) {
    if ((VAD[count] > 0.0) && (VAD[count - 2] <= 0.0)) {
      VAD[count - 2] = VAD[count] * 0.1;
      VAD[count - 1] = VAD[count] * 0.3;
      count++;
    }

    if (VAD[count] <= 0.0) {
      d = VAD[count - 1];
      if (d > 0.0) {
        VAD[count] = d * 0.3;
        VAD[count + 1] = VAD[count - 1] * 0.1;
        count += 3;
      }
    }

    count++;
  }

  end = VAD.size(1);
  for (start = 0; start < end; start++) {
    if (VAD[start] < 0.0) {
      VAD[start] = 0.0;
    }
  }

  //  fid= fopen( 'mat_vad.txt', 'wt');
  //  fprintf( fid, '%f\n', VAD);
  //  fclose( fid);
  if (LevelThresh <= 0.0) {
    LevelThresh = LevelMin;
  }

  vlen = VAD.size(1);
  logVAD.set_size(1, VAD.size(1));
  for (i = 0; i < vlen; i++) {
    logVAD[i] = 0.0;
  }

  end = VAD.size(1);
  for (start = 0; start < end; start++) {
    if (VAD[start] <= LevelThresh) {
      logVAD[start] = 0.0;
    }
  }

  c_VAD.set_size(1, VAD.size(1));
  loop_ub = VAD.size(0) * VAD.size(1);
  for (i = 0; i < loop_ub; i++) {
    c_VAD[i] = (VAD[i] > LevelThresh);
  }

  eml_find(c_VAD, r);
  VAD_greaterthan_LevelThresh.set_size(1, r.size(1));
  loop_ub = r.size(0) * r.size(1);
  for (i = 0; i < loop_ub; i++) {
    VAD_greaterthan_LevelThresh[i] = r[i];
  }

  VAD_lessthan_LevelThresh.set_size(1, VAD_greaterthan_LevelThresh.size(1));
  loop_ub = VAD_greaterthan_LevelThresh.size(0) *
    VAD_greaterthan_LevelThresh.size(1);
  for (i = 0; i < loop_ub; i++) {
    VAD_lessthan_LevelThresh[i] = VAD[static_cast<int>
      (VAD_greaterthan_LevelThresh[i]) - 1] / LevelThresh;
  }

  vlen = VAD_lessthan_LevelThresh.size(1);
  for (k = 0; k < vlen; k++) {
    VAD_lessthan_LevelThresh[k] = std::log(VAD_lessthan_LevelThresh[k]);
  }

  loop_ub = VAD_lessthan_LevelThresh.size(0) * VAD_lessthan_LevelThresh.size(1)
    - 1;
  for (i = 0; i <= loop_ub; i++) {
    logVAD[static_cast<int>(VAD_greaterthan_LevelThresh[i]) - 1] =
      VAD_lessthan_LevelThresh[i];
  }
}

static void apply_filter(const coder::array<double, 2U> &data, double
  data_Nsamples, const struct_T *glb, coder::array<double, 2U> &align_filtered)
{
  int loop_ub;
  int nx;
  double n;
  double overallGainFilter;
  double pow_of_2;
  static const double x[26] = { 0.0, 50.0, 100.0, 125.0, 160.0, 200.0, 250.0,
    300.0, 350.0, 400.0, 500.0, 600.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0,
    2000.0, 2500.0, 3000.0, 3250.0, 3500.0, 4000.0, 5000.0, 6300.0, 8000.0 };

  static const double dv[26] = { -500.0, -500.0, -500.0, -500.0, -500.0, -500.0,
    -500.0, -500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -500.0, -500.0, -500.0, -500.0, -500.0 };

  coder::array<double, 2U> b_x;
  double d;
  int loop_ub_tmp;
  coder::array<creal_T, 2U> x_fft;
  double b;
  double b_b;
  coder::array<double, 2U> y;
  coder::array<double, 2U> r;
  coder::array<double, 2U> factor;
  int j2;
  coder::array<creal_T, 2U> b_y;
  coder::array<int, 2U> r1;

  //  -------------
  // fprintf( '\tPrediction PESQ_MOS = %4.3f\n', glb.pesq_mos );
  //      clearvars -global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs glb.WHOLE_SIGNAL glb.Align_Nfft glb.Window 
  // global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs
  align_filtered.set_size(1, data.size(1));
  loop_ub = data.size(0) * data.size(1);
  for (nx = 0; nx < loop_ub; nx++) {
    align_filtered[nx] = data[nx];
  }

  n = (data_Nsamples - 150.0 * glb->Downsample) + 320.0 * (glb->Fs / 1000.0);

  //  now find the next power of 2 which is greater or equal to n
  overallGainFilter = b_log2(n);
  pow_of_2 = rt_powd_snf(2.0, std::ceil(overallGainFilter));
  overallGainFilter = interp1(x, dv);
  loop_ub = static_cast<int>(pow_of_2);
  b_x.set_size(1, loop_ub);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = 0.0;
  }

  d = 75.0 * glb->Downsample;
  loop_ub_tmp = static_cast<int>(std::floor(n - 1.0));
  for (nx = 0; nx <= loop_ub_tmp; nx++) {
    b_x[nx] = data[static_cast<int>(d + static_cast<double>(nx + 1)) - 1];
  }

  b_fft(b_x, pow_of_2, x_fft);
  b = glb->Fs / pow_of_2;
  loop_ub = static_cast<int>(pow_of_2 / 2.0 + 1.0);
  b_x.set_size(1, loop_ub);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = 0.0;
  }

  b_b = pow_of_2 / 2.0;
  if (rtIsInf(b_b) && (0.0 == b_b)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else {
    loop_ub = static_cast<int>(std::floor(b_b));
    y.set_size(1, (loop_ub + 1));
    for (nx = 0; nx <= loop_ub; nx++) {
      y[nx] = nx;
    }
  }

  nx = y.size(0) * y.size(1);
  y.set_size(1, y.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    y[nx] = y[nx] * b;
  }

  unsigned int outsize_idx_1;
  outsize_idx_1 = static_cast<unsigned int>(y.size(1));
  r.set_size(1, (static_cast<int>(outsize_idx_1)));
  loop_ub = static_cast<int>(outsize_idx_1);
  for (nx = 0; nx < loop_ub; nx++) {
    r[nx] = rtNaN;
  }

  interp1Linear(dv, y, r, x);
  loop_ub = r.size(1);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = r[nx] - overallGainFilter;
  }

  nx = b_x.size(0) * b_x.size(1);
  b_x.set_size(1, b_x.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    b_x[nx] = b_x[nx] / 20.0;
  }

  factor.set_size(1, b_x.size(1));
  nx = b_x.size(1);
  for (j2 = 0; j2 < nx; j2++) {
    factor[j2] = rt_powd_snf(10.0, b_x[j2]);
  }

  overallGainFilter = pow_of_2 / 2.0;
  if (2.0 > overallGainFilter) {
    nx = -1;
    j2 = -1;
  } else {
    nx = 0;
    j2 = static_cast<int>(overallGainFilter) - 1;
  }

  loop_ub = j2 - nx;
  b_x.set_size(1, loop_ub);
  for (j2 = 0; j2 < loop_ub; j2++) {
    b_x[j2] = factor[(nx + j2) + 1];
  }

  nx = loop_ub >> 1;
  for (int b_j1 = 0; b_j1 < nx; b_j1++) {
    j2 = (loop_ub - b_j1) - 1;
    overallGainFilter = b_x[b_j1];
    b_x[b_j1] = b_x[j2];
    b_x[j2] = overallGainFilter;
  }

  nx = factor.size(1);
  loop_ub = b_x.size(1);
  factor.set_size(factor.size(0), (factor.size(1) + b_x.size(1)));
  for (j2 = 0; j2 < loop_ub; j2++) {
    factor[nx + j2] = b_x[j2];
  }

  nx = x_fft.size(0) * x_fft.size(1);
  x_fft.set_size(1, x_fft.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    x_fft[nx].re = factor[nx] * x_fft[nx].re;
    x_fft[nx].im = factor[nx] * x_fft[nx].im;
  }

  ifft(x_fft, pow_of_2, b_y);
  if (1.0 > n) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(n);
  }

  r1.set_size(1, (loop_ub_tmp + 1));
  for (nx = 0; nx <= loop_ub_tmp; nx++) {
    r1[nx] = static_cast<int>(d + (static_cast<double>(nx) + 1.0));
  }

  for (nx = 0; nx < loop_ub; nx++) {
    align_filtered[r1[nx] - 1] = b_y[nx].re;
  }
}

static void apply_filters(const coder::array<double, 2U> &data, const struct_T
  *glb, coder::array<double, 2U> &mod_data)
{
  int sosMatrix_size_idx_0;
  int loop_ub;
  double sosMatrix_data[72];
  int i;
  int k;
  coder::array<double, 2U> b2;
  coder::array<double, 2U> a2;
  coder::array<double, 2U> A;

  //  fid= fopen( 'log_mat.txt', 'wt');
  //  fprintf( fid, '%f\n', y( 1: n));
  //  fclose( fid);
  // IIRFilt( glb.InIIR_Hsos, glb.InIIR_Nsos, data, data_Nsamples);
  //      global glb.InIIR_Hsos glb.InIIR_Nsos glb.DATAPADDING_MSEC glb.Fs
  //  now we construct the second order section matrix
  sosMatrix_size_idx_0 = glb->InIIR_Hsos.size[0];
  loop_ub = glb->InIIR_Hsos.size[0] * 6;
  if (0 <= loop_ub - 1) {
    std::memset(&sosMatrix_data[0], 0, loop_ub * sizeof(double));
  }

  //      sosMatrix_= zeros( size(glb.InIIR_Hsos,1), 6);
  //      sosMatrix( :, 4)= 1; %set a(1) to 1
  loop_ub = static_cast<int>(glb->InIIR_Nsos);
  for (i = 0; i < loop_ub; i++) {
    sosMatrix_data[i + sosMatrix_size_idx_0 * 3] = 1.0;
  }

  //  each row of sosMatrix holds [b(1*3) a(1*3)] for each section
  //       sosMatrix( :, 1: 3)= glb.InIIR_Hsos( :, 1: 3);
  //       sosMatrix( :, 5: 6)= glb.InIIR_Hsos( :, 4: 5);
  loop_ub = glb->InIIR_Hsos.size[0];
  for (i = 0; i < 3; i++) {
    for (k = 0; k < loop_ub; k++) {
      sosMatrix_data[k + sosMatrix_size_idx_0 * i] = glb->InIIR_Hsos.data[k +
        glb->InIIR_Hsos.size[0] * i];
    }
  }

  loop_ub = glb->InIIR_Hsos.size[0];
  for (i = 0; i < loop_ub; i++) {
    sosMatrix_data[i + sosMatrix_size_idx_0 * 4] = glb->InIIR_Hsos.data[i +
      glb->InIIR_Hsos.size[0] * 3];
  }

  for (i = 0; i < loop_ub; i++) {
    sosMatrix_data[i + sosMatrix_size_idx_0 * 5] = glb->InIIR_Hsos.data[i +
      glb->InIIR_Hsos.size[0] * 4];
  }

  // sosMatrix
  //  now we construct second order section direct form II filter
  //      iirdf2= dfilt.df2sos( sosMatrix);
  b2.set_size(1, 3);
  a2.set_size(1, 3);
  b2[0] = sosMatrix_data[0];
  a2[0] = sosMatrix_data[sosMatrix_size_idx_0 * 3];
  b2[1] = sosMatrix_data[sosMatrix_size_idx_0];
  a2[1] = sosMatrix_data[sosMatrix_size_idx_0 * 4];
  b2[2] = sosMatrix_data[sosMatrix_size_idx_0 * 2];
  a2[2] = sosMatrix_data[sosMatrix_size_idx_0 * 5];
  for (int m = 0; m <= sosMatrix_size_idx_0 - 2; m++) {
    int b_k;
    A.set_size(1, b2.size(1));
    loop_ub = b2.size(0) * b2.size(1);
    for (i = 0; i < loop_ub; i++) {
      A[i] = b2[i];
    }

    i = b2.size(1);
    loop_ub = b2.size(1) + 2;
    b2.set_size(1, loop_ub);
    for (k = 0; k < loop_ub; k++) {
      b2[k] = 0.0;
    }

    i--;
    for (k = 0; k < 3; k++) {
      for (b_k = 0; b_k <= i; b_k++) {
        loop_ub = k + b_k;
        b2[loop_ub] = b2[loop_ub] + sosMatrix_data[(m + sosMatrix_size_idx_0 * k)
          + 1] * A[b_k];
      }
    }

    A.set_size(1, a2.size(1));
    loop_ub = a2.size(0) * a2.size(1);
    for (i = 0; i < loop_ub; i++) {
      A[i] = a2[i];
    }

    i = a2.size(1);
    loop_ub = a2.size(1) + 2;
    a2.set_size(1, loop_ub);
    for (k = 0; k < loop_ub; k++) {
      a2[k] = 0.0;
    }

    i--;
    for (k = 0; k < 3; k++) {
      for (b_k = 0; b_k <= i; b_k++) {
        loop_ub = k + b_k;
        a2[loop_ub] = a2[loop_ub] + sosMatrix_data[(m + sosMatrix_size_idx_0 *
          (k + 3)) + 1] * A[b_k];
      }
    }
  }

  if ((b2.size(1) > 3) && (b2[b2.size(1) - 1] == 0.0)) {
    b2.set_size(b2.size(0), (b2.size(1) - 1));
  }

  if ((a2.size(1) > 3) && (a2[a2.size(1) - 1] == 0.0)) {
    a2.set_size(a2.size(0), (a2.size(1) - 1));
  }

  //      mod_data= filter( iirdf2, data);
  filter(b2, a2, data, mod_data);
}

static void apply_filters_WB(const coder::array<double, 2U> &data, const
  struct_T *glb, coder::array<double, 2U> &mod_data)
{
  double sosMatrix[6];
  coder::array<double, 2U> b_sosMatrix;
  coder::array<double, 2U> c_sosMatrix;

  //  KKW ---------
  //      global glb.WB_InIIR_Hsos glb.WB_InIIR_Nsos glb.DATAPADDING_MSEC glb.Fs 
  //  now we construct the second order section matrix
  sosMatrix[3] = 1.0;

  // set a(1) to 1
  //  each row of sosMatrix holds [b(1*3) a(1*3)] for each section
  sosMatrix[0] = glb->WB_InIIR_Hsos[0];
  sosMatrix[1] = glb->WB_InIIR_Hsos[1];
  sosMatrix[2] = glb->WB_InIIR_Hsos[2];
  sosMatrix[4] = glb->WB_InIIR_Hsos[3];
  sosMatrix[5] = glb->WB_InIIR_Hsos[4];

  // sosMatrix
  //  now we construct second order section direct form II filter
  //      iirdf2= dfilt.df2sos( sosMatrix);
  //      iirdf2_.FilterStructure = 'Direct-Form II, Second-Order Sections';
  //      iirdf2_.Arithmetic = 'double';
  //      iirdf2_.sosMatrix = sosMatrix;
  //      iirdf2_.ScaleValues = [1;1];
  //      iirdf2_.OptimizeScaleValues = true;
  //      iirdf2_.PersistentMemory = false;
  //      mod_data= filter( iirdf2, data);
  b_sosMatrix.set((&sosMatrix[0]), 1, 3);
  c_sosMatrix.set((&sosMatrix[3]), 1, 3);
  filter(b_sosMatrix, c_sosMatrix, data, mod_data);
}

static void b_apply_filter(const coder::array<double, 2U> &data, double
  data_Nsamples, const struct_T *glb, coder::array<double, 2U> &align_filtered)
{
  int loop_ub;
  int nx;
  double n;
  double overallGainFilter;
  double pow_of_2;
  static const double x[26] = { 0.0, 50.0, 100.0, 125.0, 160.0, 200.0, 250.0,
    300.0, 350.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1300.0, 1600.0,
    2000.0, 2500.0, 3000.0, 3250.0, 3500.0, 4000.0, 5000.0, 6300.0, 8000.0 };

  static const double dv[26] = { -200.0, -40.0, -20.0, -12.0, -6.0, 0.0, 4.0,
    6.0, 8.0, 10.0, 11.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0,
    12.0, 4.0, -200.0, -200.0, -200.0, -200.0 };

  coder::array<double, 2U> b_x;
  double d;
  int loop_ub_tmp;
  coder::array<creal_T, 2U> x_fft;
  double b;
  double b_b;
  coder::array<double, 2U> y;
  coder::array<double, 2U> r;
  coder::array<double, 2U> factor;
  int j2;
  coder::array<creal_T, 2U> b_y;
  coder::array<int, 2U> r1;

  //  -------------
  // fprintf( '\tPrediction PESQ_MOS = %4.3f\n', glb.pesq_mos );
  //      clearvars -global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs glb.WHOLE_SIGNAL glb.Align_Nfft glb.Window 
  // global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs
  align_filtered.set_size(1, data.size(1));
  loop_ub = data.size(0) * data.size(1);
  for (nx = 0; nx < loop_ub; nx++) {
    align_filtered[nx] = data[nx];
  }

  n = (data_Nsamples - 150.0 * glb->Downsample) + 320.0 * (glb->Fs / 1000.0);

  //  now find the next power of 2 which is greater or equal to n
  overallGainFilter = b_log2(n);
  pow_of_2 = rt_powd_snf(2.0, std::ceil(overallGainFilter));
  overallGainFilter = interp1(x, dv);
  loop_ub = static_cast<int>(pow_of_2);
  b_x.set_size(1, loop_ub);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = 0.0;
  }

  d = 75.0 * glb->Downsample;
  loop_ub_tmp = static_cast<int>(std::floor(n - 1.0));
  for (nx = 0; nx <= loop_ub_tmp; nx++) {
    b_x[nx] = data[static_cast<int>(d + static_cast<double>(nx + 1)) - 1];
  }

  b_fft(b_x, pow_of_2, x_fft);
  b = glb->Fs / pow_of_2;
  loop_ub = static_cast<int>(pow_of_2 / 2.0 + 1.0);
  b_x.set_size(1, loop_ub);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = 0.0;
  }

  b_b = pow_of_2 / 2.0;
  if (rtIsInf(b_b) && (0.0 == b_b)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else {
    loop_ub = static_cast<int>(std::floor(b_b));
    y.set_size(1, (loop_ub + 1));
    for (nx = 0; nx <= loop_ub; nx++) {
      y[nx] = nx;
    }
  }

  nx = y.size(0) * y.size(1);
  y.set_size(1, y.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    y[nx] = y[nx] * b;
  }

  unsigned int outsize_idx_1;
  outsize_idx_1 = static_cast<unsigned int>(y.size(1));
  r.set_size(1, (static_cast<int>(outsize_idx_1)));
  loop_ub = static_cast<int>(outsize_idx_1);
  for (nx = 0; nx < loop_ub; nx++) {
    r[nx] = rtNaN;
  }

  interp1Linear(dv, y, r, x);
  loop_ub = r.size(1);
  for (nx = 0; nx < loop_ub; nx++) {
    b_x[nx] = r[nx] - overallGainFilter;
  }

  nx = b_x.size(0) * b_x.size(1);
  b_x.set_size(1, b_x.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    b_x[nx] = b_x[nx] / 20.0;
  }

  factor.set_size(1, b_x.size(1));
  nx = b_x.size(1);
  for (j2 = 0; j2 < nx; j2++) {
    factor[j2] = rt_powd_snf(10.0, b_x[j2]);
  }

  overallGainFilter = pow_of_2 / 2.0;
  if (2.0 > overallGainFilter) {
    nx = -1;
    j2 = -1;
  } else {
    nx = 0;
    j2 = static_cast<int>(overallGainFilter) - 1;
  }

  loop_ub = j2 - nx;
  b_x.set_size(1, loop_ub);
  for (j2 = 0; j2 < loop_ub; j2++) {
    b_x[j2] = factor[(nx + j2) + 1];
  }

  nx = loop_ub >> 1;
  for (int b_j1 = 0; b_j1 < nx; b_j1++) {
    j2 = (loop_ub - b_j1) - 1;
    overallGainFilter = b_x[b_j1];
    b_x[b_j1] = b_x[j2];
    b_x[j2] = overallGainFilter;
  }

  nx = factor.size(1);
  loop_ub = b_x.size(1);
  factor.set_size(factor.size(0), (factor.size(1) + b_x.size(1)));
  for (j2 = 0; j2 < loop_ub; j2++) {
    factor[nx + j2] = b_x[j2];
  }

  nx = x_fft.size(0) * x_fft.size(1);
  x_fft.set_size(1, x_fft.size(1));
  loop_ub = nx - 1;
  for (nx = 0; nx <= loop_ub; nx++) {
    x_fft[nx].re = factor[nx] * x_fft[nx].re;
    x_fft[nx].im = factor[nx] * x_fft[nx].im;
  }

  ifft(x_fft, pow_of_2, b_y);
  if (1.0 > n) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(n);
  }

  r1.set_size(1, (loop_ub_tmp + 1));
  for (nx = 0; nx <= loop_ub_tmp; nx++) {
    r1[nx] = static_cast<int>(d + (static_cast<double>(nx) + 1.0));
  }

  for (nx = 0; nx < loop_ub; nx++) {
    align_filtered[r1[nx] - 1] = b_y[nx].re;
  }
}

static void b_fix_power_level(const coder::array<double, 2U> &data, double
  data_Nsamples, double maxNsamples, struct_T *glb, coder::array<double, 2U>
  &mod_data)
{
  coder::array<double, 2U> align_filtered;
  double start_point_tmp;
  double end_point_tmp;
  double end_point;
  int i;
  int vlen;
  coder::array<double, 2U> y;
  int k;

  //  done ===========
  //  this function is used for level normalization, i.e., to fix the power
  //  level of data to a preset number, and return it to mod_data.
  //      global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs
  //      global glb.TARGET_AVG_POWER
  glb->TARGET_AVG_POWER = 1.0E+7;
  apply_filter(data, data_Nsamples, glb, align_filtered);
  start_point_tmp = 75.0 * glb->Downsample;
  end_point_tmp = 320.0 * (glb->Fs / 1000.0);
  end_point = (data_Nsamples - start_point_tmp) + end_point_tmp;

  //  function pesq_testbench( testfiles, result)
  //      % used to calculate pesq score for noisy and enhanced speech
  //
  //      fid= fopen( testfiles, 'rt');
  //      fid1= fopen( result, 'wt');
  //      tline= fgetl( fid);
  //      srate= str2num( tline);
  //      % the first element is the sampling rate
  //
  //      while 1
  //          tline = fgetl(fid);
  //          if ~ischar(tline)
  //              break;
  //          end
  //          if tline== '$'
  //              % beginning of new set of clean/noisy/enhanced speech files
  //              clean= fgetl( fid); % get clean file
  //              noisy= fgetl( fid); % get noisy file
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, noisy); 
  //              noisy_pesq= pesq_measure( srate, clean, noisy);
  //              fprintf( fid1, '\nnew set of clean/noisy/enhanced speech files:\n'); 
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, noisy, noisy_pesq);
  //          elseif tline== '#'
  //              % end of testfile
  //              break;
  //          else
  //              enhanced= tline;
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, enhanced); 
  //              enhanced_pesq= pesq_measure( srate, clean, enhanced);
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, enhanced, enhanced_pesq);
  //          end
  //
  //      end
  //
  //      fclose( fid);
  //      fclose( fid1);
  //
  //
  if (start_point_tmp + 1.0 > end_point) {
    i = -1;
    vlen = 0;
  } else {
    i = static_cast<int>(start_point_tmp + 1.0) - 2;
    vlen = static_cast<int>(end_point);
  }

  vlen = (vlen - i) - 1;
  y.set_size(1, vlen);
  for (k = 0; k < vlen; k++) {
    y[k] = rt_powd_snf(align_filtered[(i + k) + 1], 2.0);
  }

  vlen = y.size(1);
  if (y.size(1) == 0) {
    start_point_tmp = 0.0;
  } else {
    start_point_tmp = y[0];
    for (k = 2; k <= vlen; k++) {
      start_point_tmp += y[k - 1];
    }
  }

  start_point_tmp = std::sqrt(1.0E+7 / (start_point_tmp / ((maxNsamples - 150.0 *
    glb->Downsample) + end_point_tmp)));

  //  fprintf( 1, '\tglobal_scale is %f\n', global_scale);
  mod_data.set_size(1, data.size(1));
  vlen = data.size(0) * data.size(1);
  for (i = 0; i < vlen; i++) {
    mod_data[i] = data[i] * start_point_tmp;
  }
}

static void compute_delay(double stop_sample, double search_range, const coder::
  array<double, 2U> &time_series1, const coder::array<double, 2U> &time_series2,
  double *best_delay, double *max_correlation)
{
  double x;
  double power_of_2;
  int i;
  coder::array<double, 2U> x1;
  int nx;
  int k;
  double h;
  double normalization;
  coder::array<double, 2U> x2;
  coder::array<double, 2U> r;
  coder::array<creal_T, 2U> y;
  coder::array<creal_T, 2U> x2_fft;
  coder::array<creal_T, 2U> b_y;

  //  fprintf( 'result_time is %f\n\n', result_time);
  x = b_log2(2.0 * ((stop_sample - 1.0) + 1.0));
  power_of_2 = rt_powd_snf(2.0, std::ceil(x));

  //  function pesq_testbench( testfiles, result)
  //      % used to calculate pesq score for noisy and enhanced speech
  //
  //      fid= fopen( testfiles, 'rt');
  //      fid1= fopen( result, 'wt');
  //      tline= fgetl( fid);
  //      srate= str2num( tline);
  //      % the first element is the sampling rate
  //
  //      while 1
  //          tline = fgetl(fid);
  //          if ~ischar(tline)
  //              break;
  //          end
  //          if tline== '$'
  //              % beginning of new set of clean/noisy/enhanced speech files
  //              clean= fgetl( fid); % get clean file
  //              noisy= fgetl( fid); % get noisy file
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, noisy); 
  //              noisy_pesq= pesq_measure( srate, clean, noisy);
  //              fprintf( fid1, '\nnew set of clean/noisy/enhanced speech files:\n'); 
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, noisy, noisy_pesq);
  //          elseif tline== '#'
  //              % end of testfile
  //              break;
  //          else
  //              enhanced= tline;
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, enhanced); 
  //              enhanced_pesq= pesq_measure( srate, clean, enhanced);
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, enhanced, enhanced_pesq);
  //          end
  //
  //      end
  //
  //      fclose( fid);
  //      fclose( fid1);
  //
  //
  //  function pesq_testbench( testfiles, result)
  //      % used to calculate pesq score for noisy and enhanced speech
  //
  //      fid= fopen( testfiles, 'rt');
  //      fid1= fopen( result, 'wt');
  //      tline= fgetl( fid);
  //      srate= str2num( tline);
  //      % the first element is the sampling rate
  //
  //      while 1
  //          tline = fgetl(fid);
  //          if ~ischar(tline)
  //              break;
  //          end
  //          if tline== '$'
  //              % beginning of new set of clean/noisy/enhanced speech files
  //              clean= fgetl( fid); % get clean file
  //              noisy= fgetl( fid); % get noisy file
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, noisy); 
  //              noisy_pesq= pesq_measure( srate, clean, noisy);
  //              fprintf( fid1, '\nnew set of clean/noisy/enhanced speech files:\n'); 
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, noisy, noisy_pesq);
  //          elseif tline== '#'
  //              % end of testfile
  //              break;
  //          else
  //              enhanced= tline;
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, enhanced); 
  //              enhanced_pesq= pesq_measure( srate, clean, enhanced);
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, enhanced, enhanced_pesq);
  //          end
  //
  //      end
  //
  //      fclose( fid);
  //      fclose( fid1);
  //
  //
  if (1.0 > stop_sample) {
    i = 0;
  } else {
    i = static_cast<int>(stop_sample);
  }

  x1.set_size(1, i);
  if (1.0 > stop_sample) {
    nx = 0;
  } else {
    nx = static_cast<int>(stop_sample);
  }

  for (k = 0; k < nx; k++) {
    x1[k] = rt_powd_snf(time_series1[k], 2.0);
  }

  nx = x1.size(1);
  if (x1.size(1) == 0) {
    x = 0.0;
  } else {
    x = x1[0];
    for (k = 2; k <= nx; k++) {
      x += x1[k - 1];
    }
  }

  if (1.0 > stop_sample) {
    i = 0;
  } else {
    i = static_cast<int>(stop_sample);
  }

  x1.set_size(1, i);
  if (1.0 > stop_sample) {
    nx = 0;
  } else {
    nx = static_cast<int>(stop_sample);
  }

  for (k = 0; k < nx; k++) {
    x1[k] = rt_powd_snf(time_series2[k], 2.0);
  }

  nx = x1.size(1);
  if (x1.size(1) == 0) {
    h = 0.0;
  } else {
    h = x1[0];
    for (k = 2; k <= nx; k++) {
      h += x1[k - 1];
    }
  }

  normalization = std::sqrt(x / ((stop_sample - 1.0) + 1.0) * ((stop_sample -
    1.0) + 1.0) / power_of_2 * (h / ((stop_sample - 1.0) + 1.0) * ((stop_sample
    - 1.0) + 1.0) / power_of_2));

  //  fprintf( 'normalization is %f\n', normalization);
  nx = static_cast<int>(power_of_2);
  x1.set_size(1, nx);
  x2.set_size(1, nx);
  for (i = 0; i < nx; i++) {
    x1[i] = 0.0;
    x2[i] = 0.0;
  }

  if (1.0 > power_of_2) {
    nx = 0;
  } else {
    nx = static_cast<int>(power_of_2);
  }

  for (i = 0; i < nx; i++) {
    x1[i] = 0.0;
  }

  if (1.0 > power_of_2) {
    nx = 0;
  } else {
    nx = static_cast<int>(power_of_2);
  }

  for (i = 0; i < nx; i++) {
    x2[i] = 0.0;
  }

  if (1.0 > stop_sample) {
    i = 0;
  } else {
    i = static_cast<int>(stop_sample);
  }

  r.set_size(1, i);
  for (k = 0; k < i; k++) {
    r[k] = std::abs(time_series1[k]);
  }

  nx = r.size(1);
  for (i = 0; i < nx; i++) {
    x1[i] = r[i];
  }

  if (1.0 > stop_sample) {
    i = 0;
  } else {
    i = static_cast<int>(stop_sample);
  }

  r.set_size(1, i);
  for (k = 0; k < i; k++) {
    r[k] = std::abs(time_series2[k]);
  }

  nx = r.size(1);
  for (i = 0; i < nx; i++) {
    x2[i] = r[i];
  }

  b_fft(x1, power_of_2, y);
  b_fft(x2, power_of_2, x2_fft);
  b_y.set_size(1, y.size(1));
  nx = y.size(0) * y.size(1);
  for (i = 0; i < nx; i++) {
    if (y[i].im == 0.0) {
      x = y[i].re / power_of_2;
      h = 0.0;
    } else if (y[i].re == 0.0) {
      x = 0.0;
      h = y[i].im / power_of_2;
    } else {
      x = y[i].re / power_of_2;
      h = y[i].im / power_of_2;
    }

    b_y[i].re = x * x2_fft[i].re - -h * x2_fft[i].im;
    b_y[i].im = x * x2_fft[i].im + -h * x2_fft[i].re;
  }

  ifft(b_y, power_of_2, y);
  *best_delay = 0.0;
  *max_correlation = 0.0;

  //  these loop can be rewritten
  i = static_cast<int>((1.0 - (-search_range)) + -1.0);
  for (k = 0; k < i; k++) {
    x = -search_range + static_cast<double>(k);
    nx = static_cast<int>((x + 1.0) + power_of_2) - 1;
    h = rt_hypotd_snf(y[nx].re, y[nx].im) / normalization;
    if (h > *max_correlation) {
      *max_correlation = h;
      *best_delay = x;
    }
  }

  i = static_cast<int>((search_range - 1.0) + 1.0);
  for (k = 0; k < i; k++) {
    h = rt_hypotd_snf(y[k].re, y[k].im) / normalization;
    if (h > *max_correlation) {
      *max_correlation = h;
      *best_delay = k;
    }
  }

  (*best_delay)--;
}

static void crude_align(const coder::array<double, 2U> &ref_logVAD, double
  ref_Nsamples, const coder::array<double, 2U> &deg_logVAD, double deg_Nsamples,
  double Utt_id, struct_T *glb)
{
  double nr;
  int nd2;
  double startr;
  double nd;
  double I_max_Y;
  double startd;
  coder::array<double, 2U> x1;
  coder::array<double, 2U> x2;
  coder::array<double, 2U> b_x;
  coder::array<creal_T, 2U> x1_fft;
  coder::array<creal_T, 2U> x2_fft;
  coder::array<creal_T, 2U> b_x1_fft;

  //      global glb.Downsample
  //      global glb.Nutterances glb.Largest_uttsize glb.Nsurf_samples glb.Crude_DelayEst 
  //      global glb.Crude_DelayConf glb.UttSearch_Start glb.UttSearch_End glb.Utt_DelayEst 
  //      global glb.Utt_Delay glb.Utt_DelayConf glb.Utt_Start glb.Utt_End
  //      global glb.MAXNUTTERANCES glb.WHOLE_SIGNAL
  //      global glb.pesq_mos glb.subj_mos glb.cond_nr
  if (Utt_id == -1.0) {
    nr = std::floor(ref_Nsamples / glb->Downsample);
    nd = std::floor(deg_Nsamples / glb->Downsample);
    startr = 1.0;
    startd = 1.0;
  } else if (Utt_id == 50.0) {
    startr = glb->UttSearch_Start[49];
    I_max_Y = glb->Utt_DelayEst[49] / glb->Downsample;
    startd = glb->UttSearch_Start[49] + I_max_Y;
    if (startd < 0.0) {
      startr = 1.0 - I_max_Y;
      startd = 1.0;
    }

    nr = glb->UttSearch_End[49] - startr;
    nd = nr;
    I_max_Y = std::floor(deg_Nsamples / glb->Downsample);
    if (startd + nr > I_max_Y) {
      nd = I_max_Y - startd;
    }

    //      fprintf( 'nr,nd is %d,%d\n', nr, nd);
  } else {
    nd2 = static_cast<int>(Utt_id) - 1;
    startr = glb->UttSearch_Start[nd2];
    I_max_Y = glb->Crude_DelayEst / glb->Downsample;
    startd = glb->UttSearch_Start[nd2] + I_max_Y;
    if (startd < 0.0) {
      startr = 1.0 - I_max_Y;
      startd = 1.0;
    }

    nr = glb->UttSearch_End[nd2] - startr;
    nd = nr;
    I_max_Y = std::floor(deg_Nsamples / glb->Downsample);
    if (startd + nr > I_max_Y + 1.0) {
      nd = (I_max_Y - startd) + 1.0;
    }
  }

  if ((1.0 > startr) || rtIsNaN(startr)) {
    startr = 1.0;
  }

  //  <- KKW
  if ((1.0 > startd) || rtIsNaN(startd)) {
    startd = 1.0;
  }

  //  <- KKW
  I_max_Y = nr;
  if ((nr > 1.0) && (nd > 1.0)) {
    double x;
    double Nx;
    int loop_ub;
    int b_j1;
    int j2;
    double ex_re;
    double ex_im;

    //  this function has other simple implementations, current implementation is 
    //  consistent with the C version
    //  % one way to do this (in time domain) =====
    //  % fprintf( 1, 'startr, nr is %d, %d\n', startr, nr);
    //  x1= ref_VAD( startr: startr+ nr- 1);
    //  x2= deg_VAD( startd: startd+ nd- 1);
    //  x1= fliplr( x1);
    //  Y= conv( x2, x1);
    //  % done =====
    //  the other way to do this (in freq domain)===
    if (!(nr > nd)) {
      I_max_Y = nd;
    }

    x = b_log2(I_max_Y);
    Nx = rt_powd_snf(2.0, std::ceil(x));
    loop_ub = static_cast<int>(2.0 * Nx);
    x1.set_size(1, loop_ub);
    x2.set_size(1, loop_ub);
    for (b_j1 = 0; b_j1 < loop_ub; b_j1++) {
      x1[b_j1] = 0.0;
      x2[b_j1] = 0.0;
    }

    if (rtIsNaN(startd)) {
      startd = 1.0;
    }

    // <<< PL: Added to avoid index 0
    if (rtIsNaN(startr)) {
      startr = 1.0;
    }

    I_max_Y = (startr + nr) - 1.0;
    if (startr > I_max_Y) {
      b_j1 = -1;
      nd2 = -1;
    } else {
      b_j1 = static_cast<int>(startr) - 2;
      nd2 = static_cast<int>(I_max_Y) - 1;
    }

    loop_ub = nd2 - b_j1;
    b_x.set_size(1, loop_ub);
    for (nd2 = 0; nd2 < loop_ub; nd2++) {
      b_x[nd2] = ref_logVAD[(b_j1 + nd2) + 1];
    }

    nd2 = loop_ub >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (loop_ub - b_j1) - 1;
      I_max_Y = b_x[b_j1];
      b_x[b_j1] = b_x[j2];
      b_x[j2] = I_max_Y;
    }

    loop_ub = b_x.size(1);
    for (b_j1 = 0; b_j1 < loop_ub; b_j1++) {
      x1[b_j1] = b_x[b_j1];
    }

    I_max_Y = (startd + nd) - 1.0;
    if (startd > I_max_Y) {
      b_j1 = 0;
      nd2 = 0;
    } else {
      b_j1 = static_cast<int>(startd) - 1;
      nd2 = static_cast<int>(I_max_Y);
    }

    loop_ub = nd2 - b_j1;
    for (nd2 = 0; nd2 < loop_ub; nd2++) {
      x2[nd2] = deg_logVAD[b_j1 + nd2];
    }

    b_fft(x1, 2.0 * Nx, x1_fft);
    b_fft(x2, 2.0 * Nx, x2_fft);
    b_x1_fft.set_size(1, x1_fft.size(1));
    loop_ub = x1_fft.size(0) * x1_fft.size(1);
    for (b_j1 = 0; b_j1 < loop_ub; b_j1++) {
      b_x1_fft[b_j1].re = x1_fft[b_j1].re * x2_fft[b_j1].re - x1_fft[b_j1].im *
        x2_fft[b_j1].im;
      b_x1_fft[b_j1].im = x1_fft[b_j1].re * x2_fft[b_j1].im + x1_fft[b_j1].im *
        x2_fft[b_j1].re;
    }

    ifft(b_x1_fft, 2.0 * Nx, x1_fft);
    j2 = (static_cast<int>(nr) + static_cast<int>(nd)) - 1;
    x2_fft.set_size(1, j2);
    for (b_j1 = 0; b_j1 < j2; b_j1++) {
      x2_fft[b_j1] = x1_fft[b_j1];
    }

    nd2 = 1;
    ex_re = x2_fft[0].re;
    ex_im = x2_fft[0].im;
    for (loop_ub = 2; loop_ub <= j2; loop_ub++) {
      boolean_T SCALEA;
      if (rtIsNaN(x2_fft[loop_ub - 1].re) || rtIsNaN(x2_fft[loop_ub - 1].im)) {
        SCALEA = false;
      } else if (rtIsNaN(ex_re) || rtIsNaN(ex_im)) {
        SCALEA = true;
      } else {
        boolean_T SCALEB;
        nd = std::abs(ex_re);
        if ((nd > 8.9884656743115785E+307) || (std::abs(ex_im) >
             8.9884656743115785E+307)) {
          SCALEA = true;
        } else {
          SCALEA = false;
        }

        startd = std::abs(x2_fft[loop_ub - 1].re);
        if ((startd > 8.9884656743115785E+307) || (std::abs(x2_fft[loop_ub - 1].
              im) > 8.9884656743115785E+307)) {
          SCALEB = true;
        } else {
          SCALEB = false;
        }

        if (SCALEA || SCALEB) {
          x = rt_hypotd_snf(ex_re / 2.0, ex_im / 2.0);
          Nx = rt_hypotd_snf(x2_fft[loop_ub - 1].re / 2.0, x2_fft[loop_ub - 1].
                             im / 2.0);
        } else {
          x = rt_hypotd_snf(ex_re, ex_im);
          Nx = rt_hypotd_snf(x2_fft[loop_ub - 1].re, x2_fft[loop_ub - 1].im);
        }

        if (iseq(x, Nx)) {
          I_max_Y = std::abs(ex_im);
          Nx = std::abs(x2_fft[loop_ub - 1].im);
          if (nd > I_max_Y) {
            startr = nd;
            nd = I_max_Y;
          } else {
            startr = I_max_Y;
          }

          if (startd > Nx) {
            I_max_Y = startd;
            startd = Nx;
          } else {
            I_max_Y = Nx;
          }

          if (startr > I_max_Y) {
            if (nd < startd) {
              x = startr - I_max_Y;
              Nx = (nd / 2.0 + startd / 2.0) / (startr / 2.0 + I_max_Y / 2.0) *
                (startd - nd);
            } else {
              x = startr;
              Nx = I_max_Y;
            }
          } else if (startr < I_max_Y) {
            if (nd > startd) {
              Nx = I_max_Y - startr;
              x = (nd / 2.0 + startd / 2.0) / (startr / 2.0 + I_max_Y / 2.0) *
                (nd - startd);
            } else {
              x = startr;
              Nx = I_max_Y;
            }
          } else {
            x = nd;
            Nx = startd;
          }

          if (iseq(x, Nx)) {
            x = rt_atan2d_snf(ex_im, ex_re);
            Nx = rt_atan2d_snf(x2_fft[loop_ub - 1].im, x2_fft[loop_ub - 1].re);
            if (iseq(x, Nx)) {
              Nx = x2_fft[loop_ub - 1].re;
              I_max_Y = x2_fft[loop_ub - 1].im;
              if (x > 0.78539816339744828) {
                if (x > 2.3561944901923448) {
                  x = -ex_im;
                  Nx = -I_max_Y;
                } else {
                  x = -ex_re;
                  Nx = -Nx;
                }
              } else if (x > -0.78539816339744828) {
                x = ex_im;
                Nx = I_max_Y;
              } else if (x > -2.3561944901923448) {
                x = ex_re;
              } else {
                x = -ex_im;
                Nx = -I_max_Y;
              }

              if (iseq(x, Nx)) {
                x = 0.0;
                Nx = 0.0;
              }
            }
          }
        }

        SCALEA = (x < Nx);
      }

      if (SCALEA) {
        ex_re = x2_fft[loop_ub - 1].re;
        ex_im = x2_fft[loop_ub - 1].im;
        nd2 = loop_ub;
      }
    }

    I_max_Y = nd2;
    if (ex_re <= 0.0) {
      I_max_Y = nr;
    }
  }

  //  fprintf( 'max_Y, I_max_Y is %f, %d\n', max_Y, I_max_Y);
  if (Utt_id == -1.0) {
    glb->Crude_DelayEst = (I_max_Y - nr) * glb->Downsample;
    glb->Crude_DelayConf = 0.0;

    //      fprintf( 1, 'I_max_Y, nr, glb.Crude_DelayEst is %f, %f, %f\n', ...
    //          I_max_Y, nr, glb.Crude_DelayEst);
  } else if (Utt_id == 50.0) {
    glb->Utt_Delay[49] = (I_max_Y - nr) * glb->Downsample + glb->Utt_DelayEst[49];

    //      fprintf( 'startr, startd, nr, nd, I_max, glb.Utt_Delay[%d] is %d, %d, %d, %d, %d, %d\n', ... 
    //            glb.MAXNUTTERANCES, startr, startd, nr, nd, ...
    //              I_max_Y, glb.Utt_Delay(glb.MAXNUTTERANCES) );
  } else {
    //      fprintf( 'I_max_Y, nr is %d, %d\n', I_max_Y, nr);
    glb->Utt_DelayEst[static_cast<int>(Utt_id) - 1] = (I_max_Y - nr) *
      glb->Downsample + glb->Crude_DelayEst;
  }
}

static void fix_power_level(const coder::array<double, 2U> &data, double
  data_Nsamples, double maxNsamples, const b_struct_T *glb, coder::array<double,
  2U> &mod_data, struct_T *b_glb)
{
  int i;
  int vlen;
  coder::array<double, 2U> align_filtered;
  double start_point_tmp;
  double end_point_tmp;
  double end_point;
  coder::array<double, 2U> y;
  int k;
  b_glb->CALIBRATE = glb->CALIBRATE;
  b_glb->Nfmax = glb->Nfmax;
  b_glb->MAXNUTTERANCES = glb->MAXNUTTERANCES;
  b_glb->MINUTTLENGTH = glb->MINUTTLENGTH;
  b_glb->WHOLE_SIGNAL = glb->WHOLE_SIGNAL;
  std::memcpy(&b_glb->UttSearch_Start[0], &glb->UttSearch_Start[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->UttSearch_End[0], &glb->UttSearch_End[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_DelayEst[0], &glb->Utt_DelayEst[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_Delay[0], &glb->Utt_Delay[0], 50U * sizeof(double));
  std::memcpy(&b_glb->Utt_DelayConf[0], &glb->Utt_DelayConf[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_Start[0], &glb->Utt_Start[0], 50U * sizeof(double));
  std::memcpy(&b_glb->Utt_End[0], &glb->Utt_End[0], 50U * sizeof(double));
  b_glb->DATAPADDING_MSEC = glb->DATAPADDING_MSEC;
  b_glb->SEARCHBUFFER = glb->SEARCHBUFFER;
  b_glb->MINSPEECHLGTH = glb->MINSPEECHLGTH;
  b_glb->JOINSPEECHLGTH = glb->JOINSPEECHLGTH;
  b_glb->Crude_DelayEst = glb->Crude_DelayEst;
  b_glb->Crude_DelayConf = glb->Crude_DelayConf;
  b_glb->Nutterances = glb->Nutterances;
  b_glb->Largest_uttsize = glb->Largest_uttsize;
  b_glb->Best_DC1 = glb->Best_DC1;
  b_glb->Best_DC2 = glb->Best_DC2;
  b_glb->Best_ED1 = glb->Best_ED1;
  b_glb->Best_D1 = glb->Best_D1;
  b_glb->Best_ED2 = glb->Best_ED2;
  b_glb->Best_D2 = glb->Best_D2;
  b_glb->Best_BP = glb->Best_BP;
  b_glb->WB_InIIR_Nsos = glb->WB_InIIR_Nsos;
  for (i = 0; i < 5; i++) {
    b_glb->WB_InIIR_Hsos[i] = glb->WB_InIIR_Hsos[i];
  }

  b_glb->Downsample = glb->Downsample;
  b_glb->Nb = glb->Nb;
  b_glb->Sl = glb->Sl;
  b_glb->Sp = glb->Sp;
  b_glb->Fs = glb->Fs;
  b_glb->InIIR_Hsos.size[0] = glb->InIIR_Hsos.size[0];
  b_glb->InIIR_Hsos.size[1] = 5;
  vlen = glb->InIIR_Hsos.size[0] * glb->InIIR_Hsos.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->InIIR_Hsos.data[0], &glb->InIIR_Hsos.data[0], vlen *
                sizeof(double));
  }

  b_glb->InIIR_Nsos = glb->InIIR_Nsos;
  b_glb->Align_Nfft = glb->Align_Nfft;
  b_glb->nr_of_hz_bands_per_bark_band.size[0] = 1;
  b_glb->nr_of_hz_bands_per_bark_band.size[1] =
    glb->nr_of_hz_bands_per_bark_band.size[1];
  vlen = glb->nr_of_hz_bands_per_bark_band.size[0] *
    glb->nr_of_hz_bands_per_bark_band.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->nr_of_hz_bands_per_bark_band.data[0],
                &glb->nr_of_hz_bands_per_bark_band.data[0], vlen * sizeof(double));
  }

  b_glb->centre_of_band_bark.size[0] = 1;
  b_glb->centre_of_band_bark.size[1] = glb->centre_of_band_bark.size[1];
  vlen = glb->centre_of_band_bark.size[0] * glb->centre_of_band_bark.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->centre_of_band_bark.data[0],
                &glb->centre_of_band_bark.data[0], vlen * sizeof(double));
  }

  b_glb->centre_of_band_hz.size[0] = 1;
  b_glb->centre_of_band_hz.size[1] = glb->centre_of_band_hz.size[1];
  vlen = glb->centre_of_band_hz.size[0] * glb->centre_of_band_hz.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->centre_of_band_hz.data[0], &glb->centre_of_band_hz.data
                [0], vlen * sizeof(double));
  }

  b_glb->width_of_band_bark.size[0] = 1;
  b_glb->width_of_band_bark.size[1] = glb->width_of_band_bark.size[1];
  vlen = glb->width_of_band_bark.size[0] * glb->width_of_band_bark.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->width_of_band_bark.data[0],
                &glb->width_of_band_bark.data[0], vlen * sizeof(double));
  }

  b_glb->width_of_band_hz.size[0] = 1;
  b_glb->width_of_band_hz.size[1] = glb->width_of_band_hz.size[1];
  vlen = glb->width_of_band_hz.size[0] * glb->width_of_band_hz.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->width_of_band_hz.data[0], &glb->width_of_band_hz.data[0],
                vlen * sizeof(double));
  }

  b_glb->pow_dens_correction_factor.size[0] = 1;
  b_glb->pow_dens_correction_factor.size[1] =
    glb->pow_dens_correction_factor.size[1];
  vlen = glb->pow_dens_correction_factor.size[0] *
    glb->pow_dens_correction_factor.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->pow_dens_correction_factor.data[0],
                &glb->pow_dens_correction_factor.data[0], vlen * sizeof(double));
  }

  b_glb->abs_thresh_power.size[0] = 1;
  b_glb->abs_thresh_power.size[1] = glb->abs_thresh_power.size[1];
  vlen = glb->abs_thresh_power.size[0] * glb->abs_thresh_power.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->abs_thresh_power.data[0], &glb->abs_thresh_power.data[0],
                vlen * sizeof(double));
  }

  b_glb->Window.size[0] = 1;
  b_glb->Window.size[1] = glb->Window.size[1];
  vlen = glb->Window.size[0] * glb->Window.size[1];
  if (0 <= vlen - 1) {
    std::memcpy(&b_glb->Window.data[0], &glb->Window.data[0], vlen * sizeof
                (double));
  }

  //  done ===========
  //  this function is used for level normalization, i.e., to fix the power
  //  level of data to a preset number, and return it to mod_data.
  //      global glb.Downsample glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs
  //      global glb.TARGET_AVG_POWER
  b_glb->TARGET_AVG_POWER = 1.0E+7;
  apply_filter(data, data_Nsamples, b_glb, align_filtered);
  start_point_tmp = 75.0 * glb->Downsample;
  end_point_tmp = 320.0 * (glb->Fs / 1000.0);
  end_point = (data_Nsamples - start_point_tmp) + end_point_tmp;

  //  function pesq_testbench( testfiles, result)
  //      % used to calculate pesq score for noisy and enhanced speech
  //
  //      fid= fopen( testfiles, 'rt');
  //      fid1= fopen( result, 'wt');
  //      tline= fgetl( fid);
  //      srate= str2num( tline);
  //      % the first element is the sampling rate
  //
  //      while 1
  //          tline = fgetl(fid);
  //          if ~ischar(tline)
  //              break;
  //          end
  //          if tline== '$'
  //              % beginning of new set of clean/noisy/enhanced speech files
  //              clean= fgetl( fid); % get clean file
  //              noisy= fgetl( fid); % get noisy file
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, noisy); 
  //              noisy_pesq= pesq_measure( srate, clean, noisy);
  //              fprintf( fid1, '\nnew set of clean/noisy/enhanced speech files:\n'); 
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, noisy, noisy_pesq);
  //          elseif tline== '#'
  //              % end of testfile
  //              break;
  //          else
  //              enhanced= tline;
  //              fprintf( 1, 'pesq_measure( %d, %s, %s)\n', srate, clean, enhanced); 
  //              enhanced_pesq= pesq_measure( srate, clean, enhanced);
  //              fprintf( fid1, '(%d, %s, %s)\t %4.3f\n', srate, ...
  //                  clean, enhanced, enhanced_pesq);
  //          end
  //
  //      end
  //
  //      fclose( fid);
  //      fclose( fid1);
  //
  //
  if (start_point_tmp + 1.0 > end_point) {
    i = -1;
    vlen = 0;
  } else {
    i = static_cast<int>(start_point_tmp + 1.0) - 2;
    vlen = static_cast<int>(end_point);
  }

  vlen = (vlen - i) - 1;
  y.set_size(1, vlen);
  for (k = 0; k < vlen; k++) {
    y[k] = rt_powd_snf(align_filtered[(i + k) + 1], 2.0);
  }

  vlen = y.size(1);
  if (y.size(1) == 0) {
    start_point_tmp = 0.0;
  } else {
    start_point_tmp = y[0];
    for (k = 2; k <= vlen; k++) {
      start_point_tmp += y[k - 1];
    }
  }

  start_point_tmp = std::sqrt(1.0E+7 / (start_point_tmp / ((maxNsamples - 150.0 *
    glb->Downsample) + end_point_tmp)));

  //  fprintf( 1, '\tglobal_scale is %f\n', global_scale);
  mod_data.set_size(1, data.size(1));
  vlen = data.size(0) * data.size(1);
  for (i = 0; i < vlen; i++) {
    mod_data[i] = data[i] * start_point_tmp;
  }
}

static void freq_resp_compensation(double number_of_frames, const coder::array<
  double, 2U> &pitch_pow_dens_ref, const double avg_pitch_pow_dens_ref_data[],
  const double avg_pitch_pow_dens_deg_data[], const c_struct_T *glb, coder::
  array<double, 2U> &mod_pitch_pow_dens_ref)
{
  int i;
  int i1;
  int loop_ub;
  int frame;

  // global glb.Nb
  i = static_cast<int>(number_of_frames);
  i1 = static_cast<int>(glb->Nb);
  mod_pitch_pow_dens_ref.set_size(i, i1);
  loop_ub = i * i1;
  for (frame = 0; frame < loop_ub; frame++) {
    mod_pitch_pow_dens_ref[frame] = 0.0;
  }

  for (loop_ub = 0; loop_ub < i1; loop_ub++) {
    double x;
    x = (avg_pitch_pow_dens_deg_data[loop_ub] + 1000.0) /
      (avg_pitch_pow_dens_ref_data[loop_ub] + 1000.0);
    if (x > 100.0) {
      x = 100.0;
    } else {
      if (x < 0.01) {
        x = 0.01;
      }
    }

    for (frame = 0; frame < i; frame++) {
      mod_pitch_pow_dens_ref[frame + mod_pitch_pow_dens_ref.size(0) * loop_ub] =
        pitch_pow_dens_ref[frame + pitch_pow_dens_ref.size(0) * loop_ub] * x;
    }
  }
}

static void freq_warping(const double hz_spectrum_data[], const c_struct_T *glb,
  double pitch_pow_dens_data[], int pitch_pow_dens_size[2])
{
  int hz_band;
  int loop_ub;
  double sum;

  //      global glb.nr_of_hz_bands_per_bark_band glb.pow_dens_correction_factor 
  //      global glb.Sp
  hz_band = 0;
  pitch_pow_dens_size[0] = 1;
  loop_ub = static_cast<int>(glb->Nb);
  pitch_pow_dens_size[1] = loop_ub;
  for (int bark_band = 0; bark_band < loop_ub; bark_band++) {
    int i;
    pitch_pow_dens_data[bark_band] = 0.0;
    sum = 0.0;
    i = static_cast<int>(glb->nr_of_hz_bands_per_bark_band.data[bark_band]);
    for (int b_i = 0; b_i < i; b_i++) {
      sum += hz_spectrum_data[hz_band];
      hz_band++;
    }

    sum *= glb->pow_dens_correction_factor.data[bark_band];
    pitch_pow_dens_data[bark_band] = sum * glb->Sp;
  }
}

static void input_filter(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &deg_data, double deg_Nsamples,
  const struct_T *glb, coder::array<double, 2U> &mod_ref_data, coder::array<
  double, 2U> &mod_deg_data)
{
  coder::array<double, 2U> b_mod_ref_data;
  coder::array<double, 2U> b_mod_deg_data;
  DC_block(ref_data, ref_Nsamples, glb, b_mod_ref_data);
  DC_block(deg_data, deg_Nsamples, glb, b_mod_deg_data);
  apply_filters(b_mod_ref_data, glb, mod_ref_data);
  apply_filters(b_mod_deg_data, glb, mod_deg_data);
}

static void intensity_warping_of(double frame, const coder::array<double, 2U>
  &pitch_pow_dens, const c_struct_T *glb, double loudness_dens_data[], int
  loudness_dens_size[2])
{
  int loop_ub;
  double h;

  // global glb.abs_thresh_power glb.Sl glb.Nb glb.centre_of_band_bark
  loudness_dens_size[0] = 1;
  loop_ub = static_cast<int>(glb->Nb);
  loudness_dens_size[1] = loop_ub;
  for (int band = 0; band < loop_ub; band++) {
    double input;
    input = pitch_pow_dens[(static_cast<int>(frame + 1.0) + pitch_pow_dens.size
      (0) * band) - 1];
    if (glb->centre_of_band_bark.data[band] < 4.0) {
      h = 6.0 / (glb->centre_of_band_bark.data[band] + 2.0);
    } else {
      h = 1.0;
    }

    if (h > 2.0) {
      h = 2.0;
    }

    h = rt_powd_snf(h, 0.15);
    h *= 0.23;
    if (input > glb->abs_thresh_power.data[band]) {
      h = rt_powd_snf(glb->abs_thresh_power.data[band] / 0.5, h) * (rt_powd_snf
        (0.5 * input / glb->abs_thresh_power.data[band] + 0.5, h) - 1.0);
    } else {
      h = 0.0;
    }

    h *= 0.1866055;
    loudness_dens_data[band] = h;
  }
}

static void multiply_with_asymmetry_factor(const double disturbance_dens_data[],
  const int disturbance_dens_size[2], double frame, const coder::array<double,
  2U> &pitch_pow_dens_ref, const coder::array<double, 2U> &pitch_pow_dens_deg,
  const c_struct_T *glb, double mod_disturbance_dens_data[], int
  mod_disturbance_dens_size[2])
{
  signed char unnamed_idx_1;
  int loop_ub;

  // global glb.Nb
  unnamed_idx_1 = static_cast<signed char>(disturbance_dens_size[1]);
  mod_disturbance_dens_size[0] = 1;
  mod_disturbance_dens_size[1] = unnamed_idx_1;
  loop_ub = unnamed_idx_1;
  if (0 <= loop_ub - 1) {
    std::memset(&mod_disturbance_dens_data[0], 0, loop_ub * sizeof(double));
  }

  loop_ub = static_cast<int>(glb->Nb);
  for (int i = 0; i < loop_ub; i++) {
    int h_tmp;
    double h;
    h_tmp = static_cast<int>(frame + 1.0) - 1;
    h = rt_powd_snf((pitch_pow_dens_deg[h_tmp + pitch_pow_dens_deg.size(0) * i]
                     + 50.0) / (pitch_pow_dens_ref[h_tmp +
      pitch_pow_dens_ref.size(0) * i] + 50.0), 1.2);
    if (h > 12.0) {
      h = 12.0;
    } else {
      if (h < 3.0) {
        h = 0.0;
      }
    }

    mod_disturbance_dens_data[i] = disturbance_dens_data[i] * h;
  }
}

static void pesq_psychoacoustic_model(const coder::array<double, 2U> &ref_data,
  double ref_Nsamples, const coder::array<double, 2U> &deg_data, double
  deg_Nsamples, const struct_T *glb, double *pesq_mos, c_struct_T *b_glb)
{
  int i;
  int max_itself_and_right;
  double maxNsamples;
  double Nf;
  double start_frame_of_bad_interval[1000];
  double stop_frame_of_bad_interval[1000];
  double start_sample_of_bad_interval[1000];
  double stop_sample_of_bad_interval[1000];
  double c_number_of_samples_in_bad_inte[1000];
  double c_delay_in_samples_in_bad_inter[1000];
  int there_is_a_bad_frame;
  double Whanning_data[514];
  int Whanning_size[1];
  double samples_to_skip_at_start;
  boolean_T exitg1;
  double d;
  double samples_to_skip_at_end;
  double y[5];
  double start_frame;
  double x_tmp_tmp;
  double x;
  int i1;
  int vlen;
  coder::array<double, 2U> ref;
  double hz_spectrum_deg_data[256];
  int loop_ub;
  coder::array<signed char, 2U> frame_is_bad;
  double sum_of_5_samples;
  coder::array<signed char, 2U> smeared_frame_is_bad;
  coder::array<double, 2U> pitch_pow_dens_ref;
  coder::array<double, 2U> pitch_pow_dens_deg;
  coder::array<double, 2U> frame_disturbance;
  coder::array<double, 2U> frame_disturbance_asym_add;
  coder::array<double, 2U> time_weight;
  coder::array<double, 2U> total_power_ref;
  int frame;
  double avg_pitch_pow_dens_ref_data[49];
  int hz_spectrum_ref_size[2];
  double avg_pitch_pow_dens_deg_data[49];
  int avg_pitch_pow_dens_deg_size[2];
  double b_Whanning_data[514];
  double hz_spectrum_ref_data[256];
  coder::array<double, 2U> b_pitch_pow_dens_ref;
  int disturbance_dens_size[2];
  int i2;
  double nn;
  double j;
  double disturbance_dens_data[49];
  int b_loop_ub;
  int c_loop_ub;
  coder::array<double, 2U> deg;
  coder::array<int, 2U> r;
  int i4;
  int tmp_data[2048];
  b_glb->CALIBRATE = glb->CALIBRATE;
  b_glb->Nfmax = glb->Nfmax;
  b_glb->MAXNUTTERANCES = glb->MAXNUTTERANCES;
  b_glb->MINUTTLENGTH = glb->MINUTTLENGTH;
  b_glb->WHOLE_SIGNAL = glb->WHOLE_SIGNAL;
  std::memcpy(&b_glb->UttSearch_Start[0], &glb->UttSearch_Start[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->UttSearch_End[0], &glb->UttSearch_End[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_DelayEst[0], &glb->Utt_DelayEst[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_Delay[0], &glb->Utt_Delay[0], 50U * sizeof(double));
  std::memcpy(&b_glb->Utt_DelayConf[0], &glb->Utt_DelayConf[0], 50U * sizeof
              (double));
  std::memcpy(&b_glb->Utt_Start[0], &glb->Utt_Start[0], 50U * sizeof(double));
  std::memcpy(&b_glb->Utt_End[0], &glb->Utt_End[0], 50U * sizeof(double));
  b_glb->DATAPADDING_MSEC = glb->DATAPADDING_MSEC;
  b_glb->SEARCHBUFFER = glb->SEARCHBUFFER;
  b_glb->MINSPEECHLGTH = glb->MINSPEECHLGTH;
  b_glb->JOINSPEECHLGTH = glb->JOINSPEECHLGTH;
  b_glb->Crude_DelayEst = glb->Crude_DelayEst;
  b_glb->Crude_DelayConf = glb->Crude_DelayConf;
  b_glb->Nutterances = glb->Nutterances;
  b_glb->Largest_uttsize = glb->Largest_uttsize;
  b_glb->Best_DC1 = glb->Best_DC1;
  b_glb->Best_DC2 = glb->Best_DC2;
  b_glb->Best_ED1 = glb->Best_ED1;
  b_glb->Best_D1 = glb->Best_D1;
  b_glb->Best_ED2 = glb->Best_ED2;
  b_glb->Best_D2 = glb->Best_D2;
  b_glb->Best_BP = glb->Best_BP;
  b_glb->WB_InIIR_Nsos = glb->WB_InIIR_Nsos;
  for (i = 0; i < 5; i++) {
    b_glb->WB_InIIR_Hsos[i] = glb->WB_InIIR_Hsos[i];
  }

  b_glb->Downsample = glb->Downsample;
  b_glb->Nb = glb->Nb;
  b_glb->Sl = glb->Sl;
  b_glb->Sp = glb->Sp;
  b_glb->Fs = glb->Fs;
  b_glb->InIIR_Hsos.size[0] = glb->InIIR_Hsos.size[0];
  b_glb->InIIR_Hsos.size[1] = 5;
  max_itself_and_right = glb->InIIR_Hsos.size[0] * glb->InIIR_Hsos.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->InIIR_Hsos.data[0], &glb->InIIR_Hsos.data[0],
                max_itself_and_right * sizeof(double));
  }

  b_glb->InIIR_Nsos = glb->InIIR_Nsos;
  b_glb->Align_Nfft = glb->Align_Nfft;
  b_glb->nr_of_hz_bands_per_bark_band.size[0] = 1;
  b_glb->nr_of_hz_bands_per_bark_band.size[1] =
    glb->nr_of_hz_bands_per_bark_band.size[1];
  max_itself_and_right = glb->nr_of_hz_bands_per_bark_band.size[0] *
    glb->nr_of_hz_bands_per_bark_band.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->nr_of_hz_bands_per_bark_band.data[0],
                &glb->nr_of_hz_bands_per_bark_band.data[0], max_itself_and_right
                * sizeof(double));
  }

  b_glb->centre_of_band_bark.size[0] = 1;
  b_glb->centre_of_band_bark.size[1] = glb->centre_of_band_bark.size[1];
  max_itself_and_right = glb->centre_of_band_bark.size[0] *
    glb->centre_of_band_bark.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->centre_of_band_bark.data[0],
                &glb->centre_of_band_bark.data[0], max_itself_and_right * sizeof
                (double));
  }

  b_glb->centre_of_band_hz.size[0] = 1;
  b_glb->centre_of_band_hz.size[1] = glb->centre_of_band_hz.size[1];
  max_itself_and_right = glb->centre_of_band_hz.size[0] *
    glb->centre_of_band_hz.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->centre_of_band_hz.data[0], &glb->centre_of_band_hz.data
                [0], max_itself_and_right * sizeof(double));
  }

  b_glb->width_of_band_bark.size[0] = 1;
  b_glb->width_of_band_bark.size[1] = glb->width_of_band_bark.size[1];
  max_itself_and_right = glb->width_of_band_bark.size[0] *
    glb->width_of_band_bark.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->width_of_band_bark.data[0],
                &glb->width_of_band_bark.data[0], max_itself_and_right * sizeof
                (double));
  }

  b_glb->width_of_band_hz.size[0] = 1;
  b_glb->width_of_band_hz.size[1] = glb->width_of_band_hz.size[1];
  max_itself_and_right = glb->width_of_band_hz.size[0] *
    glb->width_of_band_hz.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->width_of_band_hz.data[0], &glb->width_of_band_hz.data[0],
                max_itself_and_right * sizeof(double));
  }

  b_glb->pow_dens_correction_factor.size[0] = 1;
  b_glb->pow_dens_correction_factor.size[1] =
    glb->pow_dens_correction_factor.size[1];
  max_itself_and_right = glb->pow_dens_correction_factor.size[0] *
    glb->pow_dens_correction_factor.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->pow_dens_correction_factor.data[0],
                &glb->pow_dens_correction_factor.data[0], max_itself_and_right *
                sizeof(double));
  }

  b_glb->abs_thresh_power.size[0] = 1;
  b_glb->abs_thresh_power.size[1] = glb->abs_thresh_power.size[1];
  max_itself_and_right = glb->abs_thresh_power.size[0] *
    glb->abs_thresh_power.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->abs_thresh_power.data[0], &glb->abs_thresh_power.data[0],
                max_itself_and_right * sizeof(double));
  }

  b_glb->Window.size[0] = 1;
  b_glb->Window.size[1] = glb->Window.size[1];
  max_itself_and_right = glb->Window.size[0] * glb->Window.size[1];
  if (0 <= max_itself_and_right - 1) {
    std::memcpy(&b_glb->Window.data[0], &glb->Window.data[0],
                max_itself_and_right * sizeof(double));
  }

  b_glb->TARGET_AVG_POWER = glb->TARGET_AVG_POWER;

  //      global glb.CALIBRATE glb.Nfmax glb.Nb glb.Sl glb.Sp
  //      global glb.nr_of_hz_bands_per_bark_band glb.centre_of_band_bark
  //      global glb.width_of_band_hz glb.centre_of_band_hz glb.width_of_band_bark 
  //      global glb.pow_dens_correction_factor glb.abs_thresh_power
  //      global glb.Downsample glb.SEARCHBUFFER glb.DATAPADDING_MSEC glb.Fs glb.Nutterances 
  //      global glb.Utt_Start glb.Utt_End glb.Utt_Delay glb.NUMBER_OF_PSQM_FRAMES_PER_SYLLABE  
  //      global glb.Plot_Frame
  //  glb.Plot_Frame= 75; % this is the frame whose spectrum will be plotted
  b_glb->Plot_Frame = -1.0;
  b_glb->c_NUMBER_OF_PSQM_FRAMES_PER_SYL = 20.0;
  if ((ref_Nsamples > deg_Nsamples) || rtIsNaN(deg_Nsamples)) {
    maxNsamples = ref_Nsamples;
  } else {
    maxNsamples = deg_Nsamples;
  }

  Nf = glb->Downsample * 8.0;
  std::memset(&start_frame_of_bad_interval[0], 0, 1000U * sizeof(double));
  std::memset(&stop_frame_of_bad_interval[0], 0, 1000U * sizeof(double));
  std::memset(&start_sample_of_bad_interval[0], 0, 1000U * sizeof(double));
  std::memset(&stop_sample_of_bad_interval[0], 0, 1000U * sizeof(double));
  std::memset(&c_number_of_samples_in_bad_inte[0], 0, 1000U * sizeof(double));
  std::memset(&c_delay_in_samples_in_bad_inter[0], 0, 1000U * sizeof(double));
  there_is_a_bad_frame = 0;
  hann(Nf, Whanning_data, Whanning_size);
  samples_to_skip_at_start = 0.0;
  exitg1 = false;
  while ((!exitg1) && (samples_to_skip_at_start < maxNsamples / 2.0)) {
    d = samples_to_skip_at_start + 75.0 * b_glb->Downsample;
    for (max_itself_and_right = 0; max_itself_and_right < 5;
         max_itself_and_right++) {
      y[max_itself_and_right] = std::abs(ref_data[static_cast<int>(d + (
        static_cast<double>(max_itself_and_right) + 1.0)) - 1]);
    }

    if ((((y[0] + y[1]) + y[2]) + y[3]) + y[4] < 500.0) {
      samples_to_skip_at_start++;
    } else {
      exitg1 = true;
    }
  }

  //  fprintf( 'samples_to_skip_at_start is %d\n', samples_to_skip_at_start);
  samples_to_skip_at_end = 0.0;
  exitg1 = false;
  while ((!exitg1) && (samples_to_skip_at_end < maxNsamples / 2.0)) {
    d = ((maxNsamples - 75.0 * b_glb->Downsample) + 320.0 * (b_glb->Fs / 1000.0))
      - samples_to_skip_at_end;
    if (d - 4.0 > d) {
      i = 0;
      i1 = 0;
    } else {
      i = static_cast<int>(d - 4.0) - 1;
      i1 = static_cast<int>(d);
    }

    vlen = i1 - i;
    ref.set_size(1, vlen);
    for (max_itself_and_right = 0; max_itself_and_right < vlen;
         max_itself_and_right++) {
      ref[max_itself_and_right] = std::abs(ref_data[i + max_itself_and_right]);
    }

    vlen = ref.size(1);
    if (ref.size(1) == 0) {
      sum_of_5_samples = 0.0;
    } else {
      sum_of_5_samples = ref[0];
      for (max_itself_and_right = 2; max_itself_and_right <= vlen;
           max_itself_and_right++) {
        sum_of_5_samples += ref[max_itself_and_right - 1];
      }
    }

    if (sum_of_5_samples < 500.0) {
      samples_to_skip_at_end++;
    } else {
      exitg1 = true;
    }
  }

  //  fprintf( 'samples_to_skip_at_end is %d\n', samples_to_skip_at_end);
  start_frame = std::floor(samples_to_skip_at_start / (Nf / 2.0));
  x_tmp_tmp = 320.0 * (glb->Fs / 1000.0);
  x = std::floor((((maxNsamples - 150.0 * glb->Downsample) + x_tmp_tmp) -
                  samples_to_skip_at_end) / (Nf / 2.0));

  //  number of frames in speech data plus glb.DATAPADDING_MSEC
  //  fprintf( 'start/end frame is %d/%d\n', start_frame, stop_frame);
  //  fprintf( 'ref/deg power is %f/%f\n', power_ref, power_deg);
  max_itself_and_right = static_cast<int>(Nf / 2.0);
  if (0 <= max_itself_and_right - 1) {
    std::memset(&hz_spectrum_deg_data[0], 0, max_itself_and_right * sizeof
                (double));
  }

  loop_ub = static_cast<int>((x - 1.0) + 1.0);
  frame_is_bad.set_size(1, loop_ub);
  smeared_frame_is_bad.set_size(1, loop_ub);
  ref.set_size(1, loop_ub);
  for (i = 0; i < loop_ub; i++) {
    frame_is_bad[i] = 0;
    smeared_frame_is_bad[i] = 0;
    ref[i] = 0.0;
  }

  i = static_cast<int>(glb->Nb);
  pitch_pow_dens_ref.set_size(loop_ub, i);
  vlen = loop_ub * i;
  pitch_pow_dens_deg.set_size(loop_ub, i);
  for (i = 0; i < vlen; i++) {
    pitch_pow_dens_ref[i] = 0.0;
    pitch_pow_dens_deg[i] = 0.0;
  }

  frame_disturbance.set_size(1, loop_ub);
  frame_disturbance_asym_add.set_size(1, loop_ub);
  time_weight.set_size(1, loop_ub);
  total_power_ref.set_size(1, loop_ub);
  for (i = 0; i < loop_ub; i++) {
    frame_disturbance[i] = 0.0;
    frame_disturbance_asym_add[i] = 0.0;
    time_weight[i] = 0.0;
    total_power_ref[i] = 0.0;
  }

  //  fid= fopen( 'tmp_mat.txt', 'wt');
  for (frame = 0; frame < loop_ub; frame++) {
    sum_of_5_samples = (75.0 * b_glb->Downsample + 1.0) + static_cast<double>
      (frame) * (Nf / 2.0);
    vlen = Whanning_size[0];
    if (0 <= vlen - 1) {
      std::memcpy(&b_Whanning_data[0], &Whanning_data[0], vlen * sizeof(double));
    }

    short_term_fft(Nf, ref_data, b_Whanning_data, sum_of_5_samples,
                   hz_spectrum_ref_data, hz_spectrum_ref_size);
    samples_to_skip_at_start = b_glb->Nutterances;
    while ((samples_to_skip_at_start >= 1.0) && ((b_glb->Utt_Start[static_cast<
             int>(samples_to_skip_at_start) - 1] - 1.0) * b_glb->Downsample +
            1.0 > sum_of_5_samples)) {
      samples_to_skip_at_start--;
    }

    if (samples_to_skip_at_start >= 1.0) {
      samples_to_skip_at_start = b_glb->Utt_Delay[static_cast<int>
        (samples_to_skip_at_start) - 1];
    } else {
      samples_to_skip_at_start = b_glb->Utt_Delay[0];
    }

    sum_of_5_samples += samples_to_skip_at_start;
    if ((sum_of_5_samples > 0.0) && ((sum_of_5_samples + Nf) - 1.0 < maxNsamples
         + 320.0 * (b_glb->Fs / 1000.0))) {
      vlen = Whanning_size[0];
      if (0 <= vlen - 1) {
        std::memcpy(&b_Whanning_data[0], &Whanning_data[0], vlen * sizeof(double));
      }

      short_term_fft(Nf, deg_data, b_Whanning_data, sum_of_5_samples,
                     hz_spectrum_deg_data, hz_spectrum_ref_size);
    } else {
      if (0 <= max_itself_and_right - 1) {
        std::memset(&hz_spectrum_deg_data[0], 0, max_itself_and_right * sizeof
                    (double));
      }
    }

    freq_warping(hz_spectrum_ref_data, b_glb, avg_pitch_pow_dens_ref_data,
                 hz_spectrum_ref_size);
    vlen = hz_spectrum_ref_size[1];
    for (i = 0; i < vlen; i++) {
      pitch_pow_dens_ref[frame + pitch_pow_dens_ref.size(0) * i] =
        avg_pitch_pow_dens_ref_data[i];
    }

    // peak = maximum_of (pitch_pow_dens_ref, 0, glb.Nb);
    freq_warping(hz_spectrum_deg_data, b_glb, avg_pitch_pow_dens_ref_data,
                 hz_spectrum_ref_size);
    vlen = hz_spectrum_ref_size[1];
    for (i = 0; i < vlen; i++) {
      pitch_pow_dens_deg[frame + pitch_pow_dens_deg.size(0) * i] =
        avg_pitch_pow_dens_ref_data[i];
    }

    sum_of_5_samples = total_audible(static_cast<double>(frame),
      pitch_pow_dens_ref, 100.0, b_glb);
    total_audible(static_cast<double>(frame), pitch_pow_dens_deg, 100.0, b_glb);
    ref[frame] = (sum_of_5_samples < 1.0E+7);

    //      fprintf( fid, 'total_audible_pow_ref[%d] is %f\n', frame, ...
    //          total_audible_pow_ref);
  }

  //  fclose( fid);
  time_avg_audible_of((x - 1.0) + 1.0, ref, pitch_pow_dens_ref, std::floor
                      (((maxNsamples - 150.0 * glb->Downsample) + 320.0 *
                        (glb->Fs / 1000.0)) / (Nf / 2.0)) - 1.0, b_glb,
                      avg_pitch_pow_dens_ref_data, hz_spectrum_ref_size);
  time_avg_audible_of((x - 1.0) + 1.0, ref, pitch_pow_dens_deg, std::floor
                      (((maxNsamples - 150.0 * glb->Downsample) + 320.0 *
                        (glb->Fs / 1000.0)) / (Nf / 2.0)) - 1.0, b_glb,
                      avg_pitch_pow_dens_deg_data, avg_pitch_pow_dens_deg_size);

  //  fid= fopen( 'tmp_mat.txt', 'wt');
  //  fprintf( fid, '%f\n', avg_pitch_pow_dens_deg);
  //  fclose( fid);
  freq_resp_compensation((x - 1.0) + 1.0, pitch_pow_dens_ref,
    avg_pitch_pow_dens_ref_data, avg_pitch_pow_dens_deg_data, b_glb,
    b_pitch_pow_dens_ref);

  //  tmp1= pitch_pow_dens_ref';
  samples_to_skip_at_end = 1.0;
  if (0 <= loop_ub - 1) {
    disturbance_dens_size[0] = 1;
    i2 = static_cast<int>(b_glb->Nb);
  }

  for (frame = 0; frame < loop_ub; frame++) {
    sum_of_5_samples = total_audible(static_cast<double>(frame),
      b_pitch_pow_dens_ref, 1.0, b_glb);
    samples_to_skip_at_start = total_audible(static_cast<double>(frame),
      pitch_pow_dens_deg, 1.0, b_glb);
    total_power_ref[frame] = sum_of_5_samples;
    sum_of_5_samples = (sum_of_5_samples + 5000.0) / (samples_to_skip_at_start +
      5000.0);
    if (frame > 0) {
      sum_of_5_samples = 0.2 * samples_to_skip_at_end + 0.8 * sum_of_5_samples;
    }

    samples_to_skip_at_end = sum_of_5_samples;
    if (sum_of_5_samples > 5.0) {
      sum_of_5_samples = 5.0;
    } else {
      if (sum_of_5_samples < 0.0003) {
        sum_of_5_samples = 0.0003;
      }
    }

    max_itself_and_right = pitch_pow_dens_deg.size(1) - 1;
    hz_spectrum_ref_size[1] = pitch_pow_dens_deg.size(1);
    for (i = 0; i <= max_itself_and_right; i++) {
      avg_pitch_pow_dens_ref_data[i] = pitch_pow_dens_deg[frame +
        pitch_pow_dens_deg.size(0) * i] * sum_of_5_samples;
    }

    max_itself_and_right = hz_spectrum_ref_size[1];
    for (i = 0; i < max_itself_and_right; i++) {
      pitch_pow_dens_deg[frame + pitch_pow_dens_deg.size(0) * i] =
        avg_pitch_pow_dens_ref_data[i];
    }

    intensity_warping_of(static_cast<double>(frame), b_pitch_pow_dens_ref, b_glb,
                         avg_pitch_pow_dens_ref_data, hz_spectrum_ref_size);
    intensity_warping_of(static_cast<double>(frame), pitch_pow_dens_deg, b_glb,
                         avg_pitch_pow_dens_deg_data,
                         avg_pitch_pow_dens_deg_size);
    disturbance_dens_size[1] = avg_pitch_pow_dens_deg_size[1];
    max_itself_and_right = avg_pitch_pow_dens_deg_size[0] *
      avg_pitch_pow_dens_deg_size[1];
    for (i = 0; i < max_itself_and_right; i++) {
      disturbance_dens_data[i] = avg_pitch_pow_dens_deg_data[i] -
        avg_pitch_pow_dens_ref_data[i];
    }

    for (vlen = 0; vlen < i2; vlen++) {
      if ((avg_pitch_pow_dens_deg_data[vlen] < avg_pitch_pow_dens_ref_data[vlen])
          || rtIsNaN(avg_pitch_pow_dens_ref_data[vlen])) {
        sum_of_5_samples = avg_pitch_pow_dens_deg_data[vlen];
      } else {
        sum_of_5_samples = avg_pitch_pow_dens_ref_data[vlen];
      }

      d = 0.25 * sum_of_5_samples;
      if (disturbance_dens_data[vlen] > d) {
        disturbance_dens_data[vlen] -= d;

        //              disturbance_dens (band) = d- m;
      } else if (disturbance_dens_data[vlen] < -d) {
        disturbance_dens_data[vlen] += d;

        //                  disturbance_dens (band) = d+ m;
      } else {
        disturbance_dens_data[vlen] = 0.0;
      }
    }

    d = pseudo_Lp(disturbance_dens_data, 2.0, b_glb);
    frame_disturbance[frame] = d;
    if (d > 30.0) {
      there_is_a_bad_frame = 1;
    }

    multiply_with_asymmetry_factor(disturbance_dens_data, disturbance_dens_size,
      static_cast<double>(frame), b_pitch_pow_dens_ref, pitch_pow_dens_deg,
      b_glb, avg_pitch_pow_dens_ref_data, hz_spectrum_ref_size);
    frame_disturbance_asym_add[frame] = pseudo_Lp(avg_pitch_pow_dens_ref_data,
      1.0, b_glb);
  }

  //  fid= fopen( 'tmp_mat.txt', 'wt');
  //  fprintf( fid, '%f\n', frame_disturbance);
  //  fclose( fid);
  i = static_cast<int>(glb->Nutterances + -1.0);
  for (max_itself_and_right = 0; max_itself_and_right < i; max_itself_and_right
       ++) {
    sum_of_5_samples = b_glb->Utt_Delay[max_itself_and_right + 1];
    samples_to_skip_at_start = ((b_glb->Utt_Start[max_itself_and_right + 1] -
      1.0) - 75.0) * b_glb->Downsample + 1.0;
    samples_to_skip_at_end = std::floor((samples_to_skip_at_start +
      sum_of_5_samples) / (Nf / 2.0));
    j = std::floor(std::floor((((b_glb->Utt_End[max_itself_and_right] - 1.0) -
      75.0) * b_glb->Downsample + 1.0) + b_glb->Utt_Delay[max_itself_and_right])
                   / (Nf / 2.0));
    sum_of_5_samples -= b_glb->Utt_Delay[max_itself_and_right];
    if (samples_to_skip_at_end > j) {
      samples_to_skip_at_end = j;
    }

    if (samples_to_skip_at_end < 0.0) {
      samples_to_skip_at_end = 0.0;
    }

    //      fprintf( 'frame1, j, delay_jump is %d, %d, %d\n', frame1, ...
    //          j, delay_jump);
    if (sum_of_5_samples < -(Nf / 2.0)) {
      i1 = static_cast<int>((std::floor((samples_to_skip_at_start + std::abs
        (sum_of_5_samples)) / (Nf / 2.0)) + 1.0) + (1.0 - samples_to_skip_at_end));
      for (frame = 0; frame < i1; frame++) {
        nn = samples_to_skip_at_end + static_cast<double>(frame);
        if (nn < x - 1.0) {
          vlen = static_cast<int>(nn + 1.0) - 1;
          frame_disturbance[vlen] = 0.0;
          frame_disturbance_asym_add[vlen] = 0.0;
        }
      }
    }
  }

  nn = x_tmp_tmp + maxNsamples;
  max_itself_and_right = static_cast<int>(nn);
  pitch_pow_dens_ref.set_size(1, max_itself_and_right);
  for (i = 0; i < max_itself_and_right; i++) {
    pitch_pow_dens_ref[i] = 0.0;
  }

  //  fprintf( 'nn is %d\n', nn);
  d = 75.0 * glb->Downsample;
  i = static_cast<int>((nn - d) + (1.0 - (d + 1.0)));
  for (vlen = 0; vlen < i; vlen++) {
    sum_of_5_samples = (d + 1.0) + static_cast<double>(vlen);
    samples_to_skip_at_start = b_glb->Nutterances;
    while ((samples_to_skip_at_start >= 1.0) && ((b_glb->Utt_Start[static_cast<
             int>(samples_to_skip_at_start) - 1] - 1.0) * b_glb->Downsample >
            sum_of_5_samples)) {
      samples_to_skip_at_start--;
    }

    if (samples_to_skip_at_start >= 1.0) {
      samples_to_skip_at_start = b_glb->Utt_Delay[static_cast<int>
        (samples_to_skip_at_start) - 1];
    } else {
      samples_to_skip_at_start = b_glb->Utt_Delay[0];
    }

    j = sum_of_5_samples + samples_to_skip_at_start;
    samples_to_skip_at_start = 75.0 * b_glb->Downsample + 1.0;
    if (j < samples_to_skip_at_start) {
      j = samples_to_skip_at_start;
    }

    samples_to_skip_at_start = nn - 75.0 * b_glb->Downsample;
    if (j > samples_to_skip_at_start) {
      j = samples_to_skip_at_start;
    }

    pitch_pow_dens_ref[static_cast<int>(sum_of_5_samples) - 1] = deg_data[
      static_cast<int>(j) - 1];
  }

  if (there_is_a_bad_frame != 0) {
    double number_of_bad_intervals;
    for (frame = 0; frame < loop_ub; frame++) {
      frame_is_bad[frame] = static_cast<signed char>(frame_disturbance[frame] >
        30.0);
      smeared_frame_is_bad[frame] = 0;
    }

    frame_is_bad[0] = 0;
    i = static_cast<int>((((x - 1.0) - 1.0) - 2.0) + -1.0);
    for (frame = 0; frame < i; frame++) {
      signed char i3;
      i3 = frame_is_bad[frame + 2];
      vlen = i3;
      max_itself_and_right = i3;
      i1 = frame_is_bad[static_cast<int>(((static_cast<double>(frame) + 2.0) +
        1.0) + -2.0) - 1];
      if (i3 < i1) {
        vlen = i1;
      }

      if (i3 < frame_is_bad[frame + 2]) {
        max_itself_and_right = frame_is_bad[frame + 2];
      }

      i1 = frame_is_bad[static_cast<int>(((static_cast<double>(frame) + 2.0) +
        1.0) + -1.0) - 1];
      if (vlen < i1) {
        vlen = i1;
      }

      i1 = frame_is_bad[frame + 3];
      if (max_itself_and_right < i1) {
        max_itself_and_right = i1;
      }

      i1 = frame_is_bad[static_cast<int>((static_cast<double>(frame) + 2.0) +
        1.0) - 1];
      if (vlen < i1) {
        vlen = i1;
      }

      i1 = frame_is_bad[frame + 4];
      if (max_itself_and_right < i1) {
        max_itself_and_right = i1;
      }

      smeared_frame_is_bad[frame + 2] = static_cast<signed char>(vlen);
      if (vlen > max_itself_and_right) {
        smeared_frame_is_bad[frame + 2] = 0;
      }
    }

    number_of_bad_intervals = 0.0;
    nn = 0.0;
    while (nn <= x - 1.0) {
      while ((nn <= x - 1.0) && (smeared_frame_is_bad[static_cast<int>(nn + 1.0)
              - 1] == 0)) {
        nn++;
      }

      if (nn <= x - 1.0) {
        vlen = static_cast<int>(number_of_bad_intervals + 1.0) - 1;
        start_frame_of_bad_interval[vlen] = nn + 1.0;
        while ((nn <= x - 1.0) && (smeared_frame_is_bad[static_cast<int>(nn +
                 1.0) - 1] != 0)) {
          nn++;
        }

        if (nn <= x - 1.0) {
          stop_frame_of_bad_interval[vlen] = nn + 1.0;
          if (stop_frame_of_bad_interval[vlen] -
              start_frame_of_bad_interval[vlen] >= 5.0) {
            number_of_bad_intervals++;
          }
        }
      }
    }

    i = static_cast<int>((number_of_bad_intervals - 1.0) + 1.0);
    for (there_is_a_bad_frame = 0; there_is_a_bad_frame < i;
         there_is_a_bad_frame++) {
      d = 75.0 * b_glb->Downsample;
      samples_to_skip_at_start =
        ((start_frame_of_bad_interval[there_is_a_bad_frame] - 1.0) * (Nf / 2.0)
         + d) + 1.0;
      start_sample_of_bad_interval[there_is_a_bad_frame] =
        samples_to_skip_at_start;
      d += (stop_frame_of_bad_interval[there_is_a_bad_frame] - 1.0) * (Nf / 2.0)
        + Nf;
      stop_sample_of_bad_interval[there_is_a_bad_frame] = d;
      if (stop_frame_of_bad_interval[there_is_a_bad_frame] > (x - 1.0) + 1.0) {
        stop_frame_of_bad_interval[there_is_a_bad_frame] = (x - 1.0) + 1.0;
      }

      c_number_of_samples_in_bad_inte[there_is_a_bad_frame] = (d -
        samples_to_skip_at_start) + 1.0;
    }

    //      fprintf( 'number of bad intervals %d\n', number_of_bad_intervals);
    //      fprintf( '%d %d\n', number_of_samples_in_bad_interval(1), ...
    //          number_of_samples_in_bad_interval(2));
    //      fprintf( '%d %d\n', start_sample_of_bad_interval(1), ...
    //          start_sample_of_bad_interval(2));
    samples_to_skip_at_end = 4.0 * Nf;
    if (0 <= i - 1) {
      b_loop_ub = static_cast<int>(samples_to_skip_at_end);
      c_loop_ub = static_cast<int>(std::floor(samples_to_skip_at_end - 1.0));
    }

    for (there_is_a_bad_frame = 0; there_is_a_bad_frame < i;
         there_is_a_bad_frame++) {
      d = 2.0 * samples_to_skip_at_end +
        c_number_of_samples_in_bad_inte[there_is_a_bad_frame];
      max_itself_and_right = static_cast<int>(d);
      ref.set_size(1, max_itself_and_right);
      for (i1 = 0; i1 < max_itself_and_right; i1++) {
        ref[i1] = 0.0;
      }

      deg.set_size(1, max_itself_and_right);
      for (i1 = 0; i1 < max_itself_and_right; i1++) {
        deg[i1] = 0.0;
      }

      for (i1 = 0; i1 < b_loop_ub; i1++) {
        ref[i1] = 0.0;
      }

      max_itself_and_right = static_cast<int>(std::floor
        (c_number_of_samples_in_bad_inte[there_is_a_bad_frame] - 1.0));
      r.set_size(1, (max_itself_and_right + 1));
      for (i1 = 0; i1 <= max_itself_and_right; i1++) {
        r[i1] = static_cast<int>(samples_to_skip_at_end + (static_cast<double>
          (i1) + 1.0));
      }

      for (i1 = 0; i1 <= max_itself_and_right; i1++) {
        ref[r[i1] - 1] = ref_data[static_cast<int>
          (start_sample_of_bad_interval[there_is_a_bad_frame] + static_cast<
           double>(i1 + 1)) - 1];
      }

      samples_to_skip_at_start = samples_to_skip_at_end +
        c_number_of_samples_in_bad_inte[there_is_a_bad_frame];
      for (i1 = 0; i1 <= c_loop_ub; i1++) {
        tmp_data[i1] = static_cast<int>(samples_to_skip_at_start + (static_cast<
          double>(i1) + 1.0));
      }

      max_itself_and_right = static_cast<int>(std::floor(samples_to_skip_at_end
        - 1.0));
      for (i1 = 0; i1 <= max_itself_and_right; i1++) {
        ref[tmp_data[i1] - 1] = 0.0;
      }

      i1 = static_cast<int>((d - 1.0) + 1.0);
      for (vlen = 0; vlen < i1; vlen++) {
        j = (start_sample_of_bad_interval[there_is_a_bad_frame] -
             samples_to_skip_at_end) + static_cast<double>(vlen);
        nn = (maxNsamples - 75.0 * b_glb->Downsample) + 320.0 * (b_glb->Fs /
          1000.0);
        if (j <= 75.0 * b_glb->Downsample) {
          j = 75.0 * b_glb->Downsample + 1.0;
        }

        if (j > nn) {
          j = nn;
        }

        deg[vlen] = pitch_pow_dens_ref[static_cast<int>(j) - 1];
      }

      compute_delay(2.0 * samples_to_skip_at_end +
                    c_number_of_samples_in_bad_inte[there_is_a_bad_frame],
                    samples_to_skip_at_end, ref, deg,
                    &c_delay_in_samples_in_bad_inter[there_is_a_bad_frame],
                    &sum_of_5_samples);

      //          fprintf( 'delay_in_samples, best_correlation is \n\t%d, %f\n', ... 
      //              delay_in_samples, best_correlation);
      //
      if (sum_of_5_samples < 0.5) {
        c_delay_in_samples_in_bad_inter[there_is_a_bad_frame] = 0.0;
      }
    }

    if (number_of_bad_intervals > 0.0) {
      max_itself_and_right = static_cast<int>(maxNsamples + x_tmp_tmp);
      ref.set_size(1, max_itself_and_right);
      for (i1 = 0; i1 < max_itself_and_right; i1++) {
        ref[i1] = pitch_pow_dens_ref[i1];
      }

      for (there_is_a_bad_frame = 0; there_is_a_bad_frame < i;
           there_is_a_bad_frame++) {
        i1 = static_cast<int>(stop_sample_of_bad_interval[there_is_a_bad_frame]
                              + (1.0 -
          start_sample_of_bad_interval[there_is_a_bad_frame]));
        for (vlen = 0; vlen < i1; vlen++) {
          sum_of_5_samples = start_sample_of_bad_interval[there_is_a_bad_frame]
            + static_cast<double>(vlen);
          j = sum_of_5_samples +
            c_delay_in_samples_in_bad_inter[there_is_a_bad_frame];
          if (j < 1.0) {
            j = 1.0;
          }

          if (j > maxNsamples) {
            j = maxNsamples;
          }

          ref[static_cast<int>(sum_of_5_samples) - 1] = pitch_pow_dens_ref[
            static_cast<int>(j) - 1];
        }
      }

      for (there_is_a_bad_frame = 0; there_is_a_bad_frame < i;
           there_is_a_bad_frame++) {
        i1 = static_cast<int>((stop_frame_of_bad_interval[there_is_a_bad_frame]
          - 1.0) + (1.0 - start_frame_of_bad_interval[there_is_a_bad_frame]));
        for (frame = 0; frame < i1; frame++) {
          nn = (start_frame_of_bad_interval[there_is_a_bad_frame] + static_cast<
                double>(frame)) - 1.0;
          max_itself_and_right = Whanning_size[0];
          if (0 <= max_itself_and_right - 1) {
            std::memcpy(&b_Whanning_data[0], &Whanning_data[0],
                        max_itself_and_right * sizeof(double));
          }

          short_term_fft(Nf, ref, b_Whanning_data, (75.0 * b_glb->Downsample +
            nn * Nf / 2.0) + 1.0, hz_spectrum_deg_data, hz_spectrum_ref_size);
          freq_warping(hz_spectrum_deg_data, b_glb, avg_pitch_pow_dens_ref_data,
                       hz_spectrum_ref_size);
          max_itself_and_right = hz_spectrum_ref_size[1];
          for (i2 = 0; i2 < max_itself_and_right; i2++) {
            pitch_pow_dens_deg[(static_cast<int>(nn + 1.0) +
                                pitch_pow_dens_deg.size(0) * i2) - 1] =
              avg_pitch_pow_dens_ref_data[i2];
          }
        }

        samples_to_skip_at_end = 1.0;
        if (0 <= i1 - 1) {
          disturbance_dens_size[0] = 1;
          i4 = static_cast<int>(b_glb->Nb);
        }

        for (frame = 0; frame < i1; frame++) {
          nn = (start_frame_of_bad_interval[there_is_a_bad_frame] + static_cast<
                double>(frame)) - 1.0;

          //  see implementation for detail why 1 needed to be
          //  subtracted
          sum_of_5_samples = total_audible(nn, b_pitch_pow_dens_ref, 1.0, b_glb);
          samples_to_skip_at_start = total_audible(nn, pitch_pow_dens_deg, 1.0,
            b_glb);
          sum_of_5_samples = (sum_of_5_samples + 5000.0) /
            (samples_to_skip_at_start + 5000.0);
          if (nn > 0.0) {
            sum_of_5_samples = 0.2 * samples_to_skip_at_end + 0.8 *
              sum_of_5_samples;
          }

          samples_to_skip_at_end = sum_of_5_samples;
          if (sum_of_5_samples > 5.0) {
            sum_of_5_samples = 5.0;
          }

          if (sum_of_5_samples < 0.0003) {
            sum_of_5_samples = 0.0003;
          }

          max_itself_and_right = pitch_pow_dens_deg.size(1) - 1;
          hz_spectrum_ref_size[1] = pitch_pow_dens_deg.size(1);
          for (i2 = 0; i2 <= max_itself_and_right; i2++) {
            avg_pitch_pow_dens_ref_data[i2] = pitch_pow_dens_deg[(static_cast<
              int>(nn + 1.0) + pitch_pow_dens_deg.size(0) * i2) - 1] *
              sum_of_5_samples;
          }

          max_itself_and_right = hz_spectrum_ref_size[1];
          for (i2 = 0; i2 < max_itself_and_right; i2++) {
            pitch_pow_dens_deg[(static_cast<int>(nn + 1.0) +
                                pitch_pow_dens_deg.size(0) * i2) - 1] =
              avg_pitch_pow_dens_ref_data[i2];
          }

          intensity_warping_of(nn, b_pitch_pow_dens_ref, b_glb,
                               avg_pitch_pow_dens_ref_data, hz_spectrum_ref_size);
          intensity_warping_of(nn, pitch_pow_dens_deg, b_glb,
                               avg_pitch_pow_dens_deg_data,
                               avg_pitch_pow_dens_deg_size);
          disturbance_dens_size[1] = avg_pitch_pow_dens_deg_size[1];
          max_itself_and_right = avg_pitch_pow_dens_deg_size[0] *
            avg_pitch_pow_dens_deg_size[1];
          for (i2 = 0; i2 < max_itself_and_right; i2++) {
            disturbance_dens_data[i2] = avg_pitch_pow_dens_deg_data[i2] -
              avg_pitch_pow_dens_ref_data[i2];
          }

          for (vlen = 0; vlen < i4; vlen++) {
            if ((avg_pitch_pow_dens_deg_data[vlen] <
                 avg_pitch_pow_dens_ref_data[vlen]) || rtIsNaN
                (avg_pitch_pow_dens_ref_data[vlen])) {
              d = avg_pitch_pow_dens_deg_data[vlen];
            } else {
              d = avg_pitch_pow_dens_ref_data[vlen];
            }

            d *= 0.25;
            if (disturbance_dens_data[vlen] > d) {
              disturbance_dens_data[vlen] -= d;
            } else if (disturbance_dens_data[vlen] < -d) {
              disturbance_dens_data[vlen] += d;
            } else {
              disturbance_dens_data[vlen] = 0.0;
            }
          }

          sum_of_5_samples = pseudo_Lp(disturbance_dens_data, 2.0, b_glb);
          i2 = static_cast<int>(nn + 1.0) - 1;
          if ((!(frame_disturbance[i2] < sum_of_5_samples)) && (!rtIsNaN
               (sum_of_5_samples))) {
            frame_disturbance[i2] = sum_of_5_samples;
          }

          multiply_with_asymmetry_factor(disturbance_dens_data,
            disturbance_dens_size, nn, b_pitch_pow_dens_ref, pitch_pow_dens_deg,
            b_glb, avg_pitch_pow_dens_ref_data, hz_spectrum_ref_size);
          sum_of_5_samples = pseudo_Lp(avg_pitch_pow_dens_ref_data, 1.0, b_glb);
          if ((!(frame_disturbance_asym_add[i2] < sum_of_5_samples)) &&
              (!rtIsNaN(sum_of_5_samples))) {
            frame_disturbance_asym_add[i2] = sum_of_5_samples;
          }
        }
      }
    }
  }

  //  fid= fopen( 'tmp_mat1.txt', 'at');
  //  fprintf( '\n');
  for (frame = 0; frame < loop_ub; frame++) {
    time_weight[frame] = 1.0;
    if ((x - 1.0) + 1.0 > 1000.0) {
      sum_of_5_samples = std::floor((maxNsamples - 150.0 * b_glb->Downsample) /
        (Nf / 2.0));
      samples_to_skip_at_start = ((sum_of_5_samples - 1.0) - 1000.0) / 5500.0;
      if (samples_to_skip_at_start > 0.5) {
        samples_to_skip_at_start = 0.5;
      }

      time_weight[frame] = (1.0 - samples_to_skip_at_start) +
        samples_to_skip_at_start * static_cast<double>(frame) /
        (sum_of_5_samples - 1.0);
    }

    sum_of_5_samples = rt_powd_snf((total_power_ref[frame] + 100000.0) / 1.0E+7,
      0.04);

    //      if (frame== 118)
    //          fprintf( '%f\n', h);
    //          fprintf( '%f\n', frame_disturbance( 1+ frame));
    //      end
    d = frame_disturbance[frame] / sum_of_5_samples;
    frame_disturbance[frame] = d;

    //      if (frame== 118)
    //          fprintf( '%f\n', frame_disturbance( 1+ frame));
    //      end
    //
    samples_to_skip_at_start = frame_disturbance_asym_add[frame] /
      sum_of_5_samples;
    frame_disturbance_asym_add[frame] = samples_to_skip_at_start;
    if (d > 45.0) {
      frame_disturbance[frame] = 45.0;
    }

    if (samples_to_skip_at_start > 45.0) {
      frame_disturbance_asym_add[frame] = 45.0;
    }
  }

  //  fclose ( fid);
  *pesq_mos = (4.5 - 0.1 * Lpq_weight(start_frame, x - 1.0, frame_disturbance,
    time_weight)) - 0.0309 * Lpq_weight(start_frame, x - 1.0,
    frame_disturbance_asym_add, time_weight);
}

static double pseudo_Lp(const double x_data[], double p, const c_struct_T *glb)
{
  double result;
  double totalWeight;
  int i;

  //      global glb.Nb glb.width_of_band_bark
  totalWeight = 0.0;
  result = 0.0;
  i = static_cast<int>(glb->Nb + -1.0);
  for (int band = 0; band < i; band++) {
    double result_tmp;
    result_tmp = glb->width_of_band_bark.data[band + 1];
    result += rt_powd_snf(std::abs(x_data[band + 1]) * result_tmp, p);
    totalWeight += result_tmp;
  }

  return rt_powd_snf(result / totalWeight, 1.0 / p) * totalWeight;
}

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int b_u0;
    int b_u1;
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(static_cast<double>(b_u0), static_cast<double>(b_u1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

static double rt_remd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = rtNaN;
  } else if (rtIsInf(u1)) {
    y = u0;
  } else {
    double b_u1;
    if (u1 < 0.0) {
      b_u1 = std::ceil(u1);
    } else {
      b_u1 = std::floor(u1);
    }

    if ((u1 != 0.0) && (u1 != b_u1)) {
      b_u1 = std::abs(u0 / u1);
      if (!(std::abs(b_u1 - std::floor(b_u1 + 0.5)) > DBL_EPSILON * b_u1)) {
        y = 0.0 * u0;
      } else {
        y = std::fmod(u0, u1);
      }
    } else {
      y = std::fmod(u0, u1);
    }
  }

  return y;
}

static void short_term_fft(double Nf, const coder::array<double, 2U> &data,
  const double Whanning_data[], double start_sample, double hz_spectrum_data[],
  int hz_spectrum_size[2])
{
  double d;
  int y_size_idx_1;
  int k;
  int data_size[2];
  int loop_ub;
  double data_data[514];
  creal_T x1_fft_data[514];
  int x1_fft_size[2];
  coder::array<double, 2U> y;
  creal_T x_data[1024];
  double y_data[1024];
  d = (start_sample + Nf) - 1.0;
  if (start_sample > d) {
    y_size_idx_1 = 0;
    k = 0;
  } else {
    y_size_idx_1 = static_cast<int>(start_sample) - 1;
    k = static_cast<int>(d);
  }

  data_size[0] = 1;
  loop_ub = k - y_size_idx_1;
  data_size[1] = loop_ub;
  for (k = 0; k < loop_ub; k++) {
    data_data[k] = data[y_size_idx_1 + k] * Whanning_data[k];
  }

  c_fft(data_data, data_size, x1_fft_data, x1_fft_size);
  loop_ub = static_cast<int>(Nf / 2.0);
  y_size_idx_1 = static_cast<short>(loop_ub);
  for (k = 0; k < loop_ub; k++) {
    x_data[k] = x1_fft_data[k];
    y_data[k] = rt_hypotd_snf(x_data[k].re, x_data[k].im);
  }

  y.set_size(1, y_size_idx_1);
  for (k = 0; k < y_size_idx_1; k++) {
    y[k] = rt_powd_snf(y_data[k], 2.0);
  }

  hz_spectrum_size[0] = 1;
  hz_spectrum_size[1] = y.size(1);
  loop_ub = y.size(0) * y.size(1);
  for (y_size_idx_1 = 0; y_size_idx_1 < loop_ub; y_size_idx_1++) {
    hz_spectrum_data[y_size_idx_1] = y[y_size_idx_1];
  }

  hz_spectrum_data[0] = 0.0;
}

static void split_align(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &ref_logVAD, const coder::array<
  double, 2U> &deg_data, double deg_Nsamples, const coder::array<double, 2U>
  &deg_logVAD, double Utt_Start_l, double Utt_SpeechStart, double Utt_SpeechEnd,
  double Utt_End_l, double Utt_DelayEst_l, double Utt_DelayConf_l, struct_T *glb)
{
  double Utt_BPs[41];
  double Utt_ED1[41];
  double Utt_ED2[41];
  double Utt_D1[41];
  double Utt_D2[41];
  double Utt_DC1[41];
  double Utt_DC2[41];
  double Utt_Len;
  double kernel;
  double Delta;
  double Step;
  int N_BPs;
  int bp;
  int loop_ub;
  int i;
  int exitg1;
  double estdelay;
  int H_size[2];
  double H_data[1024];
  double Hsum;
  double startr;
  double startd;
  double d;
  int iindx;
  int i1;
  int count;
  coder::array<double, 2U> X1_data;
  double b_X1_data[1024];
  coder::array<creal_T, 2U> r;
  int X1_fft_size[2];
  int k;
  creal_T X1_fft_data[1024];
  coder::array<double, 2U> c_X1_data;
  coder::array<double, 2U> d_X1_data;
  creal_T X2_fft_data[1024];
  coder::array<creal_T, 2U> b_X1_fft_data;
  creal_T c_X1_fft_data[1024];
  coder::array<double, 2U> e_X1_data;
  coder::array<double, 2U> f_X1_data;
  int X1_size[2];
  coder::array<double, 2U> g_X1_data;
  coder::array<creal_T, 2U> d_X1_fft_data;
  coder::array<double, 2U> h_X1_data;
  coder::array<double, 2U> i_X1_data;
  coder::array<double, 2U> j_X1_data;
  coder::array<double, 2U> k_X1_data;
  coder::array<creal_T, 2U> e_X1_fft_data;
  coder::array<creal_T, 2U> f_X1_fft_data;
  coder::array<double, 2U> l_X1_data;
  coder::array<double, 2U> m_X1_data;

  //      global glb.MAXNUTTERANCES glb.Align_Nfft glb.Downsample glb.Window     
  //      global glb.Utt_DelayEst glb.Utt_Delay glb.UttSearch_Start glb.UttSearch_End  
  //      global glb.glb.Best_ED1 glb.Best_D1 glb.Best_DC1 glb.Best_ED2 glb.Best_D2 glb.Best_DC2 glb.Best_BP 
  std::memset(&Utt_BPs[0], 0, 41U * sizeof(double));
  std::memset(&Utt_ED1[0], 0, 41U * sizeof(double));
  std::memset(&Utt_ED2[0], 0, 41U * sizeof(double));
  std::memset(&Utt_D1[0], 0, 41U * sizeof(double));
  std::memset(&Utt_D2[0], 0, 41U * sizeof(double));
  std::memset(&Utt_DC1[0], 0, 41U * sizeof(double));
  std::memset(&Utt_DC2[0], 0, 41U * sizeof(double));
  Utt_Len = Utt_SpeechEnd - Utt_SpeechStart;
  glb->Best_DC1 = 0.0;
  glb->Best_DC2 = 0.0;
  kernel = glb->Align_Nfft / 64.0;
  Delta = glb->Align_Nfft / (4.0 * glb->Downsample);
  Step = std::floor(((0.801 * Utt_Len + 40.0 * Delta) - 1.0) / (40.0 * Delta)) *
    Delta;

  //  fprintf( 'Step is %f\n', Step);
  Delta = std::floor(Utt_Len / 10.0);
  if (Delta < 75.0) {
    Delta = 75.0;
  }

  Utt_BPs[0] = Utt_SpeechStart + Delta;
  N_BPs = -1;
  do {
    N_BPs++;
    Utt_BPs[N_BPs + 1] = Utt_BPs[N_BPs] + Step;
  } while ((Utt_BPs[N_BPs + 1] <= Utt_SpeechEnd - Delta) && (N_BPs + 2 <= 40));

  //  fprintf( 'Utt_DelayEst_l, Utt_Start_l, N_BPs is %d,%d,%d\n', ...
  //      Utt_DelayEst_l, Utt_Start_l, N_BPs);
  for (bp = 0; bp <= N_BPs; bp++) {
    glb->Utt_DelayEst[49] = Utt_DelayEst_l;
    glb->UttSearch_Start[49] = Utt_Start_l;
    glb->UttSearch_End[49] = Utt_BPs[bp];

    //      fprintf( 'bp,Utt_BPs(%d) is %d,%d\n', bp,bp,Utt_BPs(bp));
    crude_align(ref_logVAD, ref_Nsamples, deg_logVAD, deg_Nsamples, 50.0, glb);
    Utt_ED1[bp] = glb->Utt_Delay[49];
    glb->Utt_DelayEst[49] = Utt_DelayEst_l;
    glb->UttSearch_Start[49] = Utt_BPs[bp];
    glb->UttSearch_End[49] = Utt_End_l;
    crude_align(ref_logVAD, ref_Nsamples, deg_logVAD, deg_Nsamples, 50.0, glb);
    Utt_ED2[bp] = glb->Utt_Delay[49];
  }

  //  stream = fopen( 'matmat.txt', 'wt' );
  //  for count= 1: N_BPs- 1
  //      fprintf( stream, '%d\n', Utt_ED2(count));
  //  end
  //  fclose( stream );
  loop_ub = N_BPs + 1;
  for (i = 0; i < loop_ub; i++) {
    Utt_DC1[i] = -2.0;
  }

  //  stream= fopen( 'what_mmm.txt', 'at');
  do {
    exitg1 = 0;
    bp = 0;
    while ((bp + 1 <= N_BPs + 1) && (Utt_DC1[bp] > -2.0)) {
      bp++;
    }

    if (bp + 1 >= N_BPs + 2) {
      exitg1 = 1;
    } else {
      int exitg2;
      estdelay = Utt_ED1[bp];

      //      fprintf( 'bp,estdelay is %d,%d\n', bp, estdelay);
      H_size[0] = 1;
      loop_ub = static_cast<int>(glb->Align_Nfft);
      H_size[1] = loop_ub;
      if (0 <= loop_ub - 1) {
        std::memset(&H_data[0], 0, loop_ub * sizeof(double));
      }

      Hsum = 0.0;
      startr = (Utt_Start_l - 1.0) * glb->Downsample + 1.0;
      startd = startr + Utt_ED1[bp];

      //      fprintf( 'startr/startd is %d/%d\n', startr, startd);
      if (startd < 0.0) {
        startr = -Utt_ED1[bp] + 1.0;
        startd = 1.0;
      }

      if ((1.0 > startr) || rtIsNaN(startr)) {
        startr = 1.0;
      }

      //  <- KKW
      if ((1.0 > startd) || rtIsNaN(startd)) {
        startd = 1.0;
      }

      //  <- KKW
      do {
        exitg2 = 0;
        d = startd + glb->Align_Nfft;
        if (d <= deg_Nsamples + 1.0) {
          Delta = startr + glb->Align_Nfft;
          if (Delta <= (Utt_BPs[bp] - 1.0) * glb->Downsample + 1.0) {
            if (startr > Delta - 1.0) {
              i = 0;
              i1 = 0;
            } else {
              i = static_cast<int>(startr) - 1;
              i1 = static_cast<int>(Delta - 1.0);
            }

            if (startd > d - 1.0) {
              iindx = 0;
              count = 0;
            } else {
              iindx = static_cast<int>(startd) - 1;
              count = static_cast<int>(d - 1.0);
            }

            loop_ub = i1 - i;
            for (i1 = 0; i1 < loop_ub; i1++) {
              b_X1_data[i1] = ref_data[i + i1] * glb->Window.data[i1];
            }

            X1_data.set((&b_X1_data[0]), 1, loop_ub);
            b_fft(X1_data, glb->Align_Nfft, r);
            X1_fft_size[1] = r.size(1);
            loop_ub = r.size(0) * r.size(1);
            for (i = 0; i < loop_ub; i++) {
              X1_fft_data[i] = r[i];
            }

            loop_ub = count - iindx;
            for (i = 0; i < loop_ub; i++) {
              b_X1_data[i] = deg_data[iindx + i] * glb->Window.data[i];
            }

            d_X1_data.set((&b_X1_data[0]), 1, loop_ub);
            b_fft(d_X1_data, glb->Align_Nfft, r);
            loop_ub = r.size(0) * r.size(1);
            for (i = 0; i < loop_ub; i++) {
              X2_fft_data[i] = r[i];
            }

            loop_ub = X1_fft_size[1];
            for (i = 0; i < loop_ub; i++) {
              c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re -
                -X1_fft_data[i].im * X2_fft_data[i].im;
              c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im +
                -X1_fft_data[i].im * X2_fft_data[i].re;
            }

            b_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size[1]);
            ifft(b_X1_fft_data, glb->Align_Nfft, r);
            X1_fft_size[0] = 1;
            X1_fft_size[1] = r.size(1);
            loop_ub = r.size(0) * r.size(1);
            for (i = 0; i < loop_ub; i++) {
              X1_fft_data[i] = r[i];
            }

            b_abs(X1_fft_data, X1_fft_size, b_X1_data, X1_size);
            g_X1_data.set((&b_X1_data[0]), X1_size[0], X1_size[1]);
            Step = b_maximum(g_X1_data) * 0.99;
            Delta = rt_powd_snf(Step, 0.125) / kernel;

            //          fprintf( stream, '%f %f\n', v_max, n_max);
            i = static_cast<int>((glb->Align_Nfft - 1.0) + 1.0);
            for (count = 0; count < i; count++) {
              if (b_X1_data[count] > Step) {
                Hsum += Delta * kernel;
                i1 = static_cast<int>((kernel - 1.0) + (1.0 - (1.0 - kernel)));
                for (k = 0; k < i1; k++) {
                  Utt_Len = (1.0 - kernel) + static_cast<double>(k);
                  iindx = static_cast<int>(rt_remd_snf((static_cast<double>
                    (count) + Utt_Len) + glb->Align_Nfft, glb->Align_Nfft) + 1.0)
                    - 1;
                  H_data[iindx] += Delta * (kernel - std::abs(Utt_Len));
                }
              }
            }

            Delta = glb->Align_Nfft / 4.0;
            startr += Delta;
            startd += Delta;
          } else {
            exitg2 = 1;
          }
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      c_maximum(H_data, H_size, &Step, &iindx);
      Delta = iindx;
      d = glb->Align_Nfft / 2.0;
      if (static_cast<double>(iindx) - 1.0 >= d) {
        Delta = static_cast<double>(iindx) - glb->Align_Nfft;
      }

      Utt_D1[bp] = (Utt_ED1[bp] + Delta) - 1.0;
      if (Hsum > 0.0) {
        //          if (Utt_Len== 236)
        //              fprintf( 'v_max, Hsum is %f, %f\n', v_max, Hsum);
        //          end
        Utt_DC1[bp] = Step / Hsum;
      } else {
        Utt_DC1[bp] = 0.0;
      }

      //      fprintf( 'bp/startr/startd is %d/%d/%d\n', bp, startr, startd);
      while (bp + 1 < N_BPs + 1) {
        bp++;
        if ((Utt_ED1[bp] == estdelay) && (Utt_DC1[bp] <= -2.0)) {
          //              loopno= 0;
          do {
            exitg2 = 0;
            Delta = startd + glb->Align_Nfft;
            if (Delta <= deg_Nsamples + 1.0) {
              Utt_Len = startr + glb->Align_Nfft;
              if (Utt_Len <= (Utt_BPs[bp] - 1.0) * glb->Downsample + 1.0) {
                if (startr > Utt_Len - 1.0) {
                  i = 0;
                  i1 = 0;
                } else {
                  i = static_cast<int>(startr) - 1;
                  i1 = static_cast<int>(Utt_Len - 1.0);
                }

                //  %                 if (Utt_Len== 321)
                //                      fid= fopen( 'what_mat.txt', 'at');
                //                      fprintf( fid, '%f\n', glb.Window);
                //                      fclose( fid);
                //  %                     fprintf( '\n');
                //  %                 end
                if (startd > Delta - 1.0) {
                  iindx = 0;
                  count = 0;
                } else {
                  iindx = static_cast<int>(startd) - 1;
                  count = static_cast<int>(Delta - 1.0);
                }

                loop_ub = i1 - i;
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_X1_data[i1] = ref_data[i + i1] * glb->Window.data[i1];
                }

                h_X1_data.set((&b_X1_data[0]), 1, loop_ub);
                b_fft(h_X1_data, glb->Align_Nfft, r);
                X1_fft_size[1] = r.size(1);
                loop_ub = r.size(0) * r.size(1);
                for (i = 0; i < loop_ub; i++) {
                  X1_fft_data[i] = r[i];
                }

                loop_ub = count - iindx;
                for (i = 0; i < loop_ub; i++) {
                  b_X1_data[i] = deg_data[iindx + i] * glb->Window.data[i];
                }

                k_X1_data.set((&b_X1_data[0]), 1, loop_ub);
                b_fft(k_X1_data, glb->Align_Nfft, r);
                loop_ub = r.size(0) * r.size(1);
                for (i = 0; i < loop_ub; i++) {
                  X2_fft_data[i] = r[i];
                }

                loop_ub = X1_fft_size[1];
                for (i = 0; i < loop_ub; i++) {
                  c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re -
                    -X1_fft_data[i].im * X2_fft_data[i].im;
                  c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im +
                    -X1_fft_data[i].im * X2_fft_data[i].re;
                }

                f_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size[1]);
                ifft(f_X1_fft_data, glb->Align_Nfft, r);
                X1_fft_size[0] = 1;
                X1_fft_size[1] = r.size(1);
                loop_ub = r.size(0) * r.size(1);
                for (i = 0; i < loop_ub; i++) {
                  X1_fft_data[i] = r[i];
                }

                b_abs(X1_fft_data, X1_fft_size, b_X1_data, X1_size);
                m_X1_data.set((&b_X1_data[0]), X1_size[0], X1_size[1]);
                Step = 0.99 * b_maximum(m_X1_data);
                Delta = rt_powd_snf(Step, 0.125) / kernel;

                //                  fprintf( 'v_max n_max is %f %f\n', v_max, n_max); 
                i = static_cast<int>((glb->Align_Nfft - 1.0) + 1.0);
                for (count = 0; count < i; count++) {
                  if (b_X1_data[count] > Step) {
                    Hsum += Delta * kernel;
                    i1 = static_cast<int>((kernel - 1.0) + (1.0 - (1.0 - kernel)));
                    for (k = 0; k < i1; k++) {
                      Utt_Len = (1.0 - kernel) + static_cast<double>(k);
                      iindx = static_cast<int>(rt_remd_snf((static_cast<double>
                        (count) + Utt_Len) + glb->Align_Nfft, glb->Align_Nfft) +
                        1.0) - 1;
                      H_data[iindx] += Delta * (kernel - std::abs(Utt_Len));
                    }
                  }
                }

                Delta = glb->Align_Nfft / 4.0;
                startr += Delta;
                startd += Delta;

                //                  loopno= loopno+ 1;
              } else {
                exitg2 = 1;
              }
            } else {
              exitg2 = 1;
            }
          } while (exitg2 == 0);

          //              fprintf( 'loopno is %d\n', loopno);
          c_maximum(H_data, H_size, &Step, &iindx);
          Delta = iindx;

          //              fprintf( 'I_max is %d ', I_max);
          if (static_cast<double>(iindx) - 1.0 >= d) {
            Delta = static_cast<double>(iindx) - glb->Align_Nfft;
          }

          Utt_D1[bp] = (estdelay + Delta) - 1.0;
          if (Hsum > 0.0) {
            //                  fprintf( 'v_max Hsum is %f %f\n', v_max, Hsum);
            Utt_DC1[bp] = Step / Hsum;
          } else {
            Utt_DC1[bp] = 0.0;
          }
        }
      }
    }
  } while (exitg1 == 0);

  //  fclose( stream);
  for (bp = 0; bp <= N_BPs; bp++) {
    if (Utt_DC1[bp] > Utt_DelayConf_l) {
      Utt_DC2[bp] = -2.0;
    } else {
      Utt_DC2[bp] = 0.0;
    }
  }

  H_size[0] = 1;
  loop_ub = static_cast<int>(glb->Align_Nfft);
  H_size[1] = loop_ub;
  if (0 <= loop_ub - 1) {
    std::memset(&H_data[0], 0, loop_ub * sizeof(double));
  }

  do {
    exitg1 = 0;
    bp = N_BPs;
    while ((bp + 1 >= 1) && (Utt_DC2[bp] > -2.0)) {
      bp--;
    }

    if (bp + 1 < 1) {
      exitg1 = 1;
    } else {
      estdelay = Utt_ED2[bp];
      if (0 <= loop_ub - 1) {
        std::memset(&H_data[0], 0, loop_ub * sizeof(double));
      }

      Hsum = 0.0;
      startr = ((Utt_End_l - 1.0) * glb->Downsample + 1.0) - glb->Align_Nfft;
      startd = startr + Utt_ED2[bp];

      //      fprintf( '***NEW startr is %d\n', startr);
      //      fprintf( 'startr/d, deg_Nsamples is %d/%d, %d\n', startr,startd, ... 
      //          deg_Nsamples);
      //      fprintf( 'deg_data has %d elements\n', numel( deg_data));
      if (startd + glb->Align_Nfft > deg_Nsamples + 1.0) {
        startd = (deg_Nsamples - glb->Align_Nfft) + 1.0;
        startr = startd - Utt_ED2[bp];
      }

      while ((startd >= 1.0) && (startr >= (Utt_BPs[bp] - 1.0) * glb->Downsample
              + 1.0)) {
        d = (startr + glb->Align_Nfft) - 1.0;
        if (startr > d) {
          i = 0;
          i1 = 0;
        } else {
          i = static_cast<int>(startr) - 1;
          i1 = static_cast<int>(d);
        }

        d = (startd + glb->Align_Nfft) - 1.0;
        if (startd > d) {
          iindx = 0;
          count = 0;
        } else {
          iindx = static_cast<int>(startd) - 1;
          count = static_cast<int>(d);
        }

        k = i1 - i;
        for (i1 = 0; i1 < k; i1++) {
          b_X1_data[i1] = ref_data[i + i1] * glb->Window.data[i1];
        }

        c_X1_data.set((&b_X1_data[0]), 1, k);
        b_fft(c_X1_data, glb->Align_Nfft, r);
        X1_fft_size[1] = r.size(1);
        k = r.size(0) * r.size(1);
        for (i = 0; i < k; i++) {
          X1_fft_data[i] = r[i];
        }

        k = count - iindx;
        for (i = 0; i < k; i++) {
          b_X1_data[i] = deg_data[iindx + i] * glb->Window.data[i];
        }

        e_X1_data.set((&b_X1_data[0]), 1, k);
        b_fft(e_X1_data, glb->Align_Nfft, r);
        k = r.size(0) * r.size(1);
        for (i = 0; i < k; i++) {
          X2_fft_data[i] = r[i];
        }

        k = X1_fft_size[1];
        for (i = 0; i < k; i++) {
          c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re -
            -X1_fft_data[i].im * X2_fft_data[i].im;
          c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im +
            -X1_fft_data[i].im * X2_fft_data[i].re;
        }

        d_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size[1]);
        ifft(d_X1_fft_data, glb->Align_Nfft, r);
        X1_fft_size[0] = 1;
        X1_fft_size[1] = r.size(1);
        k = r.size(0) * r.size(1);
        for (i = 0; i < k; i++) {
          X1_fft_data[i] = r[i];
        }

        b_abs(X1_fft_data, X1_fft_size, b_X1_data, X1_size);
        j_X1_data.set((&b_X1_data[0]), X1_size[0], X1_size[1]);
        Step = b_maximum(j_X1_data) * 0.99;
        Delta = rt_powd_snf(Step, 0.125) / kernel;
        i = static_cast<int>((glb->Align_Nfft - 1.0) + 1.0);
        for (count = 0; count < i; count++) {
          if (b_X1_data[count] > Step) {
            Hsum += Delta * kernel;
            i1 = static_cast<int>((kernel - 1.0) + (1.0 - (1.0 - kernel)));
            for (k = 0; k < i1; k++) {
              Utt_Len = (1.0 - kernel) + static_cast<double>(k);
              iindx = static_cast<int>(rt_remd_snf((static_cast<double>(count) +
                Utt_Len) + glb->Align_Nfft, glb->Align_Nfft) + 1.0) - 1;
              H_data[iindx] += Delta * (kernel - std::abs(Utt_Len));
            }
          }
        }

        Delta = glb->Align_Nfft / 4.0;
        startr -= Delta;
        startd -= Delta;
      }

      c_maximum(H_data, H_size, &Step, &iindx);
      Delta = iindx;
      d = glb->Align_Nfft / 2.0;
      if (static_cast<double>(iindx) - 1.0 >= d) {
        Delta = static_cast<double>(iindx) - glb->Align_Nfft;
      }

      Utt_D2[bp] = (Utt_ED2[bp] + Delta) - 1.0;
      if (Hsum > 0.0) {
        Utt_DC2[bp] = Step / Hsum;
      } else {
        Utt_DC2[bp] = 0.0;
      }

      while (bp + 1 > 1) {
        bp--;
        if ((Utt_ED2[bp] == estdelay) && (Utt_DC2[bp] <= -2.0)) {
          while ((startd >= 1.0) && (startr >= (Utt_BPs[bp] - 1.0) *
                  glb->Downsample + 1.0)) {
            Delta = (startr + glb->Align_Nfft) - 1.0;
            if (startr > Delta) {
              i = 0;
              i1 = 0;
            } else {
              i = static_cast<int>(startr) - 1;
              i1 = static_cast<int>(Delta);
            }

            Delta = (startd + glb->Align_Nfft) - 1.0;
            if (startd > Delta) {
              iindx = 0;
              count = 0;
            } else {
              iindx = static_cast<int>(startd) - 1;
              count = static_cast<int>(Delta);
            }

            k = i1 - i;
            for (i1 = 0; i1 < k; i1++) {
              b_X1_data[i1] = ref_data[i + i1] * glb->Window.data[i1];
            }

            f_X1_data.set((&b_X1_data[0]), 1, k);
            b_fft(f_X1_data, glb->Align_Nfft, r);
            X1_fft_size[1] = r.size(1);
            k = r.size(0) * r.size(1);
            for (i = 0; i < k; i++) {
              X1_fft_data[i].re = r[i].re;
              X1_fft_data[i].im = -r[i].im;
            }

            k = count - iindx;
            for (i = 0; i < k; i++) {
              b_X1_data[i] = deg_data[iindx + i] * glb->Window.data[i];
            }

            i_X1_data.set((&b_X1_data[0]), 1, k);
            b_fft(i_X1_data, glb->Align_Nfft, r);
            k = r.size(0) * r.size(1);
            for (i = 0; i < k; i++) {
              X2_fft_data[i] = r[i];
            }

            k = X1_fft_size[1];
            for (i = 0; i < k; i++) {
              c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re -
                X1_fft_data[i].im * X2_fft_data[i].im;
              c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im +
                X1_fft_data[i].im * X2_fft_data[i].re;
            }

            e_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size[1]);
            ifft(e_X1_fft_data, glb->Align_Nfft, r);
            X1_fft_size[0] = 1;
            X1_fft_size[1] = r.size(1);
            k = r.size(0) * r.size(1);
            for (i = 0; i < k; i++) {
              X1_fft_data[i] = r[i];
            }

            b_abs(X1_fft_data, X1_fft_size, b_X1_data, X1_size);
            l_X1_data.set((&b_X1_data[0]), X1_size[0], X1_size[1]);
            Step = b_maximum(l_X1_data) * 0.99;
            Delta = rt_powd_snf(Step, 0.125) / kernel;
            i = static_cast<int>((glb->Align_Nfft - 1.0) + 1.0);
            for (count = 0; count < i; count++) {
              if (b_X1_data[count] > Step) {
                Hsum += Delta * kernel;
                i1 = static_cast<int>((kernel - 1.0) + (1.0 - (1.0 - kernel)));
                for (k = 0; k < i1; k++) {
                  Utt_Len = (1.0 - kernel) + static_cast<double>(k);
                  iindx = static_cast<int>(rt_remd_snf((static_cast<double>
                    (count) + Utt_Len) + glb->Align_Nfft, glb->Align_Nfft) + 1.0)
                    - 1;
                  H_data[iindx] += Delta * (kernel - std::abs(Utt_Len));
                }
              }
            }

            Delta = glb->Align_Nfft / 4.0;
            startr -= Delta;
            startd -= Delta;
          }

          c_maximum(H_data, H_size, &Step, &iindx);
          Delta = iindx;
          if (static_cast<double>(iindx) - 1.0 >= d) {
            Delta = static_cast<double>(iindx) - glb->Align_Nfft;
          }

          Utt_D2[bp] = (estdelay + Delta) - 1.0;
          if (Hsum > 0.0) {
            Utt_DC2[bp] = Step / Hsum;
          } else {
            Utt_DC2[bp] = 0.0;
          }
        }
      }
    }
  } while (exitg1 == 0);

  //  fid= fopen( 'uttinfo_mat.txt', 'wt');
  //  fprintf( fid, '%f\n', Utt_D2);
  //  fprintf( fid, '\n');
  //  fprintf( fid, '%f\n', Utt_DC2);
  //  fclose( fid);
  //  fprintf( 'Utt_Len, N_BPs is %d, %d\n', Utt_Len, N_BPs);
  for (bp = 0; bp <= N_BPs; bp++) {
    if ((std::abs(Utt_D2[bp] - Utt_D1[bp]) >= glb->Downsample) && (Utt_DC1[bp] +
         Utt_DC2[bp] > glb->Best_DC1 + glb->Best_DC2) && (Utt_DC1[bp] >
         Utt_DelayConf_l) && (Utt_DC2[bp] > Utt_DelayConf_l)) {
      glb->Best_ED1 = Utt_ED1[bp];
      glb->Best_D1 = Utt_D1[bp];
      glb->Best_DC1 = Utt_DC1[bp];
      glb->Best_ED2 = Utt_ED2[bp];
      glb->Best_D2 = Utt_D2[bp];
      glb->Best_DC2 = Utt_DC2[bp];
      glb->Best_BP = Utt_BPs[bp];

      //          fprintf( 'in loop...');
    }
  }
}

static void time_align(const coder::array<double, 2U> &ref_data, const coder::
  array<double, 2U> &deg_data, double deg_Nsamples, double Utt_id, struct_T *glb)
{
  int loop_ub;
  double H_data[1024];
  int startr_tmp;
  double startr;
  double startd;
  double d;
  double Hsum;
  int i;
  int X1_fft_size_idx_1;
  int k;
  double kernel;
  int b_loop_ub;
  int nx;
  double X2_data[1024];
  coder::array<double, 2U> b_H_data;
  coder::array<double, 2U> b_X2_data;
  coder::array<creal_T, 2U> r;
  coder::array<double, 2U> c_X2_data;
  creal_T X1_fft_data[1024];
  coder::array<double, 2U> d_X2_data;
  creal_T X2_fft_data[1024];
  coder::array<creal_T, 2U> b_X1_fft_data;
  creal_T c_X1_fft_data[1024];
  coder::array<creal_T, 2U> d_X1_fft_data;
  boolean_T exitg2;
  coder::array<boolean_T, 2U> x;
  coder::array<short, 2U> ii;
  double c_H_data[1024];

  //  if (Utt_Len== 236)
  //      fid= fopen( 'matmat.txt', 'wt');
  //      fprintf( fid, 'N_BPs is %d\n', N_BPs);
  //      fprintf( fid, 'glb.Utt_DelayConf is %f\n', Utt_DelayConf_l);
  //      fprintf( fid, 'ED2\t ED1\t D2\t D1\t DC2\t DC1\t BPs\n');
  //      for bp= 1: N_BPs- 1
  //          fprintf( fid, '%d\t %d\t %d\t %d\t %f\t %f\t %d\n', Utt_ED2( bp), ... 
  //              Utt_ED1( bp), Utt_D2(bp), Utt_D1(bp), Utt_DC2(bp),...
  //              Utt_DC1( bp), Utt_BPs( bp));
  //      end
  //      fclose( fid);
  //  end
  //      global glb.Utt_DelayEst glb.Utt_Delay glb.Utt_DelayConf glb.UttSearch_Start glb.UttSearch_End  
  //      global glb.Align_Nfft glb.Downsample glb.Window
  loop_ub = static_cast<int>(glb->Align_Nfft);
  if (0 <= loop_ub - 1) {
    std::memset(&H_data[0], 0, loop_ub * sizeof(double));
  }

  startr_tmp = static_cast<int>(Utt_id) - 1;
  startr = (glb->UttSearch_Start[startr_tmp] - 1.0) * glb->Downsample + 1.0;
  startd = startr + glb->Utt_DelayEst[startr_tmp];
  if (startd < 0.0) {
    startr = 1.0 - glb->Utt_DelayEst[startr_tmp];
    startd = 1.0;
  }

  int exitg1;
  do {
    exitg1 = 0;
    d = startd + glb->Align_Nfft;
    if (d <= deg_Nsamples) {
      Hsum = startr + glb->Align_Nfft;
      if (Hsum <= (glb->UttSearch_End[startr_tmp] - 1.0) * glb->Downsample) {
        if (startr > Hsum - 1.0) {
          i = 0;
          X1_fft_size_idx_1 = 0;
        } else {
          i = static_cast<int>(startr) - 1;
          X1_fft_size_idx_1 = static_cast<int>(Hsum - 1.0);
        }

        if (startd > d - 1.0) {
          nx = 0;
          k = 0;
        } else {
          nx = static_cast<int>(startd) - 1;
          k = static_cast<int>(d - 1.0);
        }

        //  find cross-correlation between X1 and X2
        b_loop_ub = X1_fft_size_idx_1 - i;
        for (X1_fft_size_idx_1 = 0; X1_fft_size_idx_1 < b_loop_ub;
             X1_fft_size_idx_1++) {
          X2_data[X1_fft_size_idx_1] = ref_data[i + X1_fft_size_idx_1] *
            glb->Window.data[X1_fft_size_idx_1];
        }

        b_X2_data.set((&X2_data[0]), 1, b_loop_ub);
        b_fft(b_X2_data, glb->Align_Nfft, r);
        X1_fft_size_idx_1 = r.size(1);
        b_loop_ub = r.size(0) * r.size(1);
        for (i = 0; i < b_loop_ub; i++) {
          X1_fft_data[i] = r[i];
        }

        b_loop_ub = k - nx;
        for (i = 0; i < b_loop_ub; i++) {
          X2_data[i] = deg_data[nx + i] * glb->Window.data[i];
        }

        d_X2_data.set((&X2_data[0]), 1, b_loop_ub);
        b_fft(d_X2_data, glb->Align_Nfft, r);
        b_loop_ub = r.size(0) * r.size(1);
        for (i = 0; i < b_loop_ub; i++) {
          X2_fft_data[i] = r[i];
        }

        for (i = 0; i < X1_fft_size_idx_1; i++) {
          c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re -
            -X1_fft_data[i].im * X2_fft_data[i].im;
          c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im +
            -X1_fft_data[i].im * X2_fft_data[i].re;
        }

        d_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size_idx_1);
        ifft(d_X1_fft_data, glb->Align_Nfft, r);
        X1_fft_size_idx_1 = r.size(1);
        b_loop_ub = r.size(0) * r.size(1);
        for (i = 0; i < b_loop_ub; i++) {
          X1_fft_data[i] = r[i];
        }

        nx = static_cast<short>(X1_fft_size_idx_1);
        for (k = 0; k < X1_fft_size_idx_1; k++) {
          X2_data[k] = rt_hypotd_snf(X1_fft_data[k].re, X1_fft_data[k].im);
        }

        if (nx <= 2) {
          if ((X2_data[0] < X2_data[1]) || (rtIsNaN(X2_data[0]) && (!rtIsNaN
                (X2_data[1])))) {
            kernel = X2_data[1];
          } else {
            kernel = X2_data[0];
          }
        } else {
          if (!rtIsNaN(X2_data[0])) {
            b_loop_ub = 1;
          } else {
            b_loop_ub = 0;
            k = 2;
            exitg2 = false;
            while ((!exitg2) && (k <= nx)) {
              if (!rtIsNaN(X2_data[k - 1])) {
                b_loop_ub = k;
                exitg2 = true;
              } else {
                k++;
              }
            }
          }

          if (b_loop_ub == 0) {
            kernel = X2_data[0];
          } else {
            kernel = X2_data[b_loop_ub - 1];
            i = b_loop_ub + 1;
            for (k = i; k <= nx; k++) {
              if (kernel < X2_data[k - 1]) {
                kernel = X2_data[k - 1];
              }
            }
          }
        }

        Hsum = kernel * 0.99;
        x.set_size(1, nx);
        for (i = 0; i < nx; i++) {
          x[i] = (X2_data[i] > Hsum);
        }

        nx = x.size(1);
        b_loop_ub = 0;
        ii.set_size(1, x.size(1));
        k = 0;
        exitg2 = false;
        while ((!exitg2) && (k <= nx - 1)) {
          if (x[k]) {
            b_loop_ub++;
            ii[b_loop_ub - 1] = static_cast<short>(k + 1);
            if (b_loop_ub >= nx) {
              exitg2 = true;
            } else {
              k++;
            }
          } else {
            k++;
          }
        }

        if (1 > b_loop_ub) {
          b_loop_ub = 0;
        }

        ii.set_size(ii.size(0), b_loop_ub);
        nx = ii.size(1);
        b_loop_ub = ii.size(0) * ii.size(1);
        for (i = 0; i < b_loop_ub; i++) {
          X2_data[i] = ii[i];
        }

        Hsum = rt_powd_snf(Hsum, 0.125);
        for (i = 0; i < nx; i++) {
          c_H_data[i] = H_data[static_cast<int>(X2_data[i]) - 1] + Hsum;
        }

        for (i = 0; i < nx; i++) {
          H_data[static_cast<int>(X2_data[i]) - 1] = c_H_data[i];
        }

        Hsum = glb->Align_Nfft / 4.0;
        startr += Hsum;
        startd += Hsum;
      } else {
        exitg1 = 1;
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (loop_ub == 0) {
    Hsum = 0.0;
  } else {
    Hsum = H_data[0];
    for (k = 2; k <= loop_ub; k++) {
      Hsum += H_data[k - 1];
    }
  }

  kernel = glb->Align_Nfft / 64.0;
  b_loop_ub = static_cast<int>(kernel);
  if (0 <= b_loop_ub - 1) {
    std::memset(&X2_data[0], 0, b_loop_ub * sizeof(double));
  }

  X2_data[0] = 1.0;
  i = static_cast<int>(kernel + -1.0);
  for (nx = 0; nx < i; nx++) {
    d = 1.0 - ((static_cast<double>(nx) + 2.0) - 1.0) / kernel;
    X2_data[nx + 1] = d;
    X2_data[static_cast<int>((glb->Align_Nfft - (static_cast<double>(nx) + 2.0))
      + 2.0) - 1] = d;
  }

  b_H_data.set((&H_data[0]), 1, loop_ub);
  b_fft(b_H_data, glb->Align_Nfft, r);
  X1_fft_size_idx_1 = r.size(1);
  loop_ub = r.size(0) * r.size(1);
  for (i = 0; i < loop_ub; i++) {
    X1_fft_data[i] = r[i];
  }

  c_X2_data.set((&X2_data[0]), 1, b_loop_ub);
  b_fft(c_X2_data, glb->Align_Nfft, r);
  loop_ub = r.size(0) * r.size(1);
  for (i = 0; i < loop_ub; i++) {
    X2_fft_data[i] = r[i];
  }

  for (i = 0; i < X1_fft_size_idx_1; i++) {
    c_X1_fft_data[i].re = X1_fft_data[i].re * X2_fft_data[i].re - X1_fft_data[i]
      .im * X2_fft_data[i].im;
    c_X1_fft_data[i].im = X1_fft_data[i].re * X2_fft_data[i].im + X1_fft_data[i]
      .im * X2_fft_data[i].re;
  }

  b_X1_fft_data.set((&c_X1_fft_data[0]), 1, X1_fft_size_idx_1);
  ifft(b_X1_fft_data, glb->Align_Nfft, r);
  X1_fft_size_idx_1 = r.size(1);
  loop_ub = r.size(0) * r.size(1);
  for (i = 0; i < loop_ub; i++) {
    X1_fft_data[i] = r[i];
  }

  if (Hsum > 0.0) {
    nx = static_cast<short>(X1_fft_size_idx_1);
    for (k = 0; k < X1_fft_size_idx_1; k++) {
      H_data[k] = rt_hypotd_snf(X1_fft_data[k].re, X1_fft_data[k].im);
    }

    loop_ub = nx - 1;
    for (i = 0; i <= loop_ub; i++) {
      H_data[i] /= Hsum;
    }
  } else {
    nx = 1;
    H_data[0] = 0.0;
  }

  if (nx <= 2) {
    if (nx == 1) {
      kernel = H_data[0];
      b_loop_ub = 1;
    } else if ((H_data[0] < H_data[1]) || (rtIsNaN(H_data[0]) && (!rtIsNaN
                 (H_data[1])))) {
      kernel = H_data[1];
      b_loop_ub = 2;
    } else {
      kernel = H_data[0];
      b_loop_ub = 1;
    }
  } else {
    if (!rtIsNaN(H_data[0])) {
      b_loop_ub = 1;
    } else {
      b_loop_ub = 0;
      k = 2;
      exitg2 = false;
      while ((!exitg2) && (k <= nx)) {
        if (!rtIsNaN(H_data[k - 1])) {
          b_loop_ub = k;
          exitg2 = true;
        } else {
          k++;
        }
      }
    }

    if (b_loop_ub == 0) {
      kernel = H_data[0];
      b_loop_ub = 1;
    } else {
      kernel = H_data[b_loop_ub - 1];
      i = b_loop_ub + 1;
      for (k = i; k <= nx; k++) {
        if (kernel < H_data[k - 1]) {
          kernel = H_data[k - 1];
          b_loop_ub = k;
        }
      }
    }
  }

  Hsum = b_loop_ub;
  if (static_cast<double>(b_loop_ub) - 1.0 >= glb->Align_Nfft / 2.0) {
    Hsum = static_cast<double>(b_loop_ub) - glb->Align_Nfft;
  }

  glb->Utt_Delay[startr_tmp] = (glb->Utt_DelayEst[startr_tmp] + Hsum) - 1.0;
  glb->Utt_DelayConf[startr_tmp] = kernel;
}

static void time_avg_audible_of(double number_of_frames, const coder::array<
  double, 2U> &silent, const coder::array<double, 2U> &pitch_pow_dens, double
  total_number_of_frames, const c_struct_T *glb, double avg_pitch_pow_dens_data[],
  int avg_pitch_pow_dens_size[2])
{
  int loop_ub;
  double result;

  //      global glb.Nb glb.abs_thresh_power
  avg_pitch_pow_dens_size[0] = 1;
  loop_ub = static_cast<int>(glb->Nb);
  avg_pitch_pow_dens_size[1] = loop_ub;
  for (int band = 0; band < loop_ub; band++) {
    int i;
    avg_pitch_pow_dens_data[band] = 0.0;
    result = 0.0;
    i = static_cast<int>(number_of_frames);
    for (int frame = 0; frame < i; frame++) {
      if (!(silent[frame] != 0.0)) {
        double h;
        h = pitch_pow_dens[frame + pitch_pow_dens.size(0) * band];
        if (h > 100.0 * glb->abs_thresh_power.data[band]) {
          result += h;
        }
      }

      avg_pitch_pow_dens_data[band] = result / total_number_of_frames;
    }
  }
}

static double total_audible(double frame, const coder::array<double, 2U>
  &pitch_pow_dens, double factor, const c_struct_T *glb)
{
  double total_audible_pow;
  int i;

  //      global glb.Nb glb.abs_thresh_power
  total_audible_pow = 0.0;
  i = static_cast<int>(glb->Nb + -1.0);
  for (int band = 0; band < i; band++) {
    double h;
    h = pitch_pow_dens[(static_cast<int>(frame + 1.0) + pitch_pow_dens.size(0) *
                        (band + 1)) - 1];
    if (h > factor * glb->abs_thresh_power.data[band + 1]) {
      total_audible_pow += h;
    }
  }

  return total_audible_pow;
}

static void utterance_locate(const coder::array<double, 2U> &ref_data, double
  ref_Nsamples, const coder::array<double, 2U> &ref_VAD, const coder::array<
  double, 2U> &ref_logVAD, const coder::array<double, 2U> &deg_data, double
  deg_Nsamples, const coder::array<double, 2U> &deg_logVAD, struct_T *glb)
{
  int Utt_num;
  int speech_flag;
  double this_start;
  double last_end;
  double Utt_SpeechEnd;
  int i;
  int count;
  double Utt_Start_l;
  int Utt_id;
  coder::array<double, 2U> x;
  double varargin_1[50];
  boolean_T exitg1;

  //  confidence
  //      global glb.Nutterances glb.Utt_Delay glb.Utt_DelayConf glb.Utt_Start glb.Utt_End glb.Utt_DelayEst 
  //      global glb.MINUTTLENGTH glb.Downsample glb.SEARCHBUFFER
  //      global glb.Crude_DelayEst glb.Nutterances glb.UttSearch_Start glb.UttSearch_End 
  Utt_num = 0;
  speech_flag = 0;
  this_start = std::floor(ref_Nsamples / glb->Downsample);
  last_end = 50.0 - glb->Crude_DelayEst / glb->Downsample;
  Utt_SpeechEnd = std::floor((deg_Nsamples - glb->Crude_DelayEst) /
    glb->Downsample);
  i = static_cast<int>(this_start);
  for (count = 0; count < i; count++) {
    Utt_Start_l = ref_VAD[count];
    if ((Utt_Start_l > 0.0) && (speech_flag == 0)) {
      speech_flag = 1;
      glb->UttSearch_Start[Utt_num] = (static_cast<double>(count) + 1.0) - 75.0;

      //          if( glb.UttSearch_Start(Utt_num)< 0 )
      //              glb.UttSearch_Start(Utt_num)= 0;
      //          end
      if (glb->UttSearch_Start[Utt_num] < 1.0) {
        glb->UttSearch_Start[Utt_num] = 1.0;
      }
    }

    if (((Utt_Start_l == 0.0) || (static_cast<double>(count) + 1.0 == this_start
          - 1.0)) && (speech_flag == 1)) {
      speech_flag = 0;
      glb->UttSearch_End[Utt_num] = (static_cast<double>(count) + 1.0) + 75.0;

      //          if( glb.UttSearch_End(Utt_num) > VAD_length - 1 )
      //              glb.UttSearch_End(Utt_num) = VAD_length -1;
      //          end
      if (glb->UttSearch_End[Utt_num] > this_start) {
        glb->UttSearch_End[Utt_num] = this_start;
      }

      if (((static_cast<double>(count) + 1.0) - (static_cast<double>(count) +
            1.0) >= 50.0) && (static_cast<double>(count) + 1.0 < Utt_SpeechEnd -
           50.0) && (static_cast<double>(count) + 1.0 > last_end)) {
        Utt_num++;
      }
    }
  }

  glb->Nutterances = static_cast<double>(Utt_num + 1) - 1.0;
  for (Utt_id = 0; Utt_id < Utt_num; Utt_id++) {
    // fprintf( 1, 'Utt_id is %d\n', Utt_id);
    crude_align(ref_logVAD, ref_Nsamples, deg_logVAD, deg_Nsamples, static_cast<
                double>(Utt_id) + 1.0, glb);
    time_align(ref_data, deg_data, deg_Nsamples, static_cast<double>(Utt_id) +
               1.0, glb);
  }

  //  fprintf( 1, 'glb.Nutterances is %d\n', glb.Nutterances);
  //  fid= fopen( 'mat_utt.txt', 'wt');
  //  fprintf( fid, '%d\n', glb.UttSearch_Start( 1: glb.Nutterances));
  //  fprintf( fid, '\n');
  //  fprintf( fid, '%d\n', glb.UttSearch_End( 1: glb.Nutterances));
  //  fclose(fid);
  //      global glb.Largest_uttsize glb.MINUTTLENGTH glb.Crude_DelayEst
  //      global glb.Downsample glb.SEARCHBUFFER glb.Nutterances glb.Utt_Start
  //      global glb.Utt_End glb.Utt_Delay
  Utt_num = 0;
  speech_flag = 0;
  Utt_Start_l = std::floor(ref_Nsamples / glb->Downsample);

  //  fprintf( 1, 'VAD_length is %d\n', VAD_length);
  this_start = 50.0 - glb->Crude_DelayEst / glb->Downsample;
  last_end = std::floor((deg_Nsamples - glb->Crude_DelayEst) / glb->Downsample);
  i = static_cast<int>(Utt_Start_l);
  for (count = 0; count < i; count++) {
    Utt_SpeechEnd = ref_VAD[count];
    if ((Utt_SpeechEnd > 0.0) && (speech_flag == 0)) {
      speech_flag = 1;
      glb->Utt_Start[Utt_num] = static_cast<double>(count) + 1.0;
    }

    if (((Utt_SpeechEnd == 0.0) || (static_cast<double>(count) + 1.0 ==
          Utt_Start_l)) && (speech_flag == 1)) {
      speech_flag = 0;
      glb->Utt_End[Utt_num] = static_cast<double>(count) + 1.0;
      if (((static_cast<double>(count) + 1.0) - (static_cast<double>(count) +
            1.0) >= 50.0) && (static_cast<double>(count) + 1.0 < last_end - 50.0)
          && (static_cast<double>(count) + 1.0 > this_start)) {
        Utt_num++;
      }
    }
  }

  glb->Utt_Start[0] = 76.0;
  if ((1.0 > glb->Nutterances) || rtIsNaN(glb->Nutterances)) {
    glb->Nutterances = 1.0;
  }

  // <<< PL: Added to avoid index 0
  glb->Utt_End[static_cast<int>(glb->Nutterances) - 1] = (Utt_Start_l - 75.0) +
    1.0;
  i = static_cast<int>(glb->Nutterances + -1.0);
  for (Utt_num = 0; Utt_num < i; Utt_num++) {
    this_start = std::floor(((glb->Utt_Start[Utt_num + 1] - 1.0) + (glb->
      Utt_End[Utt_num] - 1.0)) / 2.0);
    glb->Utt_Start[Utt_num + 1] = this_start + 1.0;
    glb->Utt_End[Utt_num] = this_start + 1.0;
  }

  if ((glb->Utt_Start[0] - 1.0) * glb->Downsample + glb->Utt_Delay[0] < 75.0 *
      glb->Downsample) {
    glb->Utt_Start[0] = (std::floor(((glb->Downsample - 1.0) - glb->Utt_Delay[0])
      / glb->Downsample) + 75.0) + 1.0;
  }

  //  fprintf( 'glb.Utt_End(%d) is %d\n', glb.Nutterances, glb.Utt_End(glb.Nutterances)); 
  //  fprintf( 'last_end is %d\n', last_end);
  //  fprintf( 'glb.Utt_Delay(%d) is %d\n', glb.Nutterances, glb.Utt_Delay(glb.Nutterances)); 
  if (((glb->Utt_End[static_cast<int>(glb->Nutterances) - 1] - 1.0) *
       glb->Downsample + 1.0) + glb->Utt_Delay[static_cast<int>(glb->Nutterances)
      - 1] > (deg_Nsamples - 75.0 * glb->Downsample) + 1.0) {
    glb->Utt_End[static_cast<int>(glb->Nutterances) - 1] = (std::floor
      ((deg_Nsamples - glb->Utt_Delay[static_cast<int>(glb->Nutterances) - 1]) /
       glb->Downsample) - 75.0) + 1.0;
  }

  i = static_cast<int>(glb->Nutterances + -1.0);
  for (Utt_num = 0; Utt_num < i; Utt_num++) {
    Utt_SpeechEnd = glb->Utt_Delay[Utt_num + 1];
    this_start = (glb->Utt_Start[Utt_num + 1] - 1.0) * glb->Downsample +
      Utt_SpeechEnd;
    last_end = (glb->Utt_End[Utt_num] - 1.0) * glb->Downsample + glb->
      Utt_Delay[Utt_num];
    if (this_start < last_end) {
      this_start = std::floor((this_start + last_end) / 2.0);
      glb->Utt_Start[Utt_num + 1] = std::floor((((glb->Downsample - 1.0) +
        this_start) - Utt_SpeechEnd) / glb->Downsample) + 1.0;
      glb->Utt_End[Utt_num] = std::floor((this_start - glb->Utt_Delay[Utt_num]) /
        glb->Downsample) + 1.0;
    }
  }

  x.set_size(1, 50);
  for (i = 0; i < 50; i++) {
    this_start = glb->Utt_End[i] - glb->Utt_Start[i];
    varargin_1[i] = this_start;
    x[i] = this_start;
  }

  if (!rtIsNaN(x[0])) {
    speech_flag = 1;
  } else {
    speech_flag = 0;
    count = 2;
    exitg1 = false;
    while ((!exitg1) && (count <= 50)) {
      if (!rtIsNaN(x[count - 1])) {
        speech_flag = count;
        exitg1 = true;
      } else {
        count++;
      }
    }
  }

  if (speech_flag == 0) {
    glb->Largest_uttsize = varargin_1[0];
  } else {
    this_start = varargin_1[speech_flag - 1];
    i = speech_flag + 1;
    for (count = i; count < 51; count++) {
      Utt_Start_l = varargin_1[count - 1];
      if (this_start < Utt_Start_l) {
        this_start = Utt_Start_l;
      }
    }

    glb->Largest_uttsize = this_start;
  }

  //  fid= fopen( 'mat_utt_info.txt', 'wt');
  //  fprintf( fid, 'glb.Utt_DelayEst: \n');
  //  fprintf( fid, '%d\n', glb.Utt_DelayEst( 1: glb.Nutterances));
  //  fprintf( fid, 'glb.Utt_Delay:\n');
  //  fprintf( fid, '%d\n', glb.Utt_Delay(1: glb.Nutterances));
  //  fprintf( fid, 'glb.Utt_Delay confidence:\n');
  //  fprintf( fid, '%f\n', glb.Utt_DelayConf(1: glb.Nutterances));
  //  fprintf( fid, 'glb.Utt_Start: \n');
  //  fprintf( fid, '%d\n', glb.Utt_Start( 1: glb.Nutterances));
  //  fprintf( fid, 'glb.Utt_End: \n');
  //  fprintf( fid, '%d\n', glb.Utt_End(1: glb.Nutterances));
  //  fclose( fid);
  //      global glb.Nutterances glb.MAXNUTTERANCES glb.Downsample glb.SEARCHBUFFER 
  //      global glb.Utt_DelayEst glb.Utt_Delay glb.Utt_DelayConf glb.UttSearch_Start 
  //      global glb.Utt_Start glb.Utt_End glb.Largest_uttsize glb.UttSearch_End 
  //      global glb.Best_ED1 glb.Best_D1 glb.Best_DC1 glb.Best_ED2 glb.Best_D2 glb.Best_DC2 glb.Best_BP 
  Utt_id = 0;
  int exitg2;
  do {
    exitg2 = 0;
    i = Utt_id + 1;
    if ((i <= glb->Nutterances) && (glb->Nutterances <= 50.0)) {
      double Utt_End_l;
      this_start = glb->Utt_DelayConf[Utt_id];
      Utt_Start_l = glb->Utt_Start[Utt_id];
      Utt_End_l = glb->Utt_End[Utt_id];
      last_end = glb->Utt_Start[Utt_id];
      if (!(1.0 < last_end)) {
        last_end = 1.0;
      }

      //  <- KKW
      // fprintf( 'SpeechStart is %d\n', Utt_SpeechStart);
      while ((last_end < Utt_End_l) && (ref_VAD[static_cast<int>(last_end) - 1] <=
              0.0)) {
        last_end++;
      }

      // find the SpeechStart for each utterance
      Utt_SpeechEnd = glb->Utt_End[Utt_id];

      //      fprintf( 'SpeechEnd is %d\n', Utt_SpeechEnd);
      while ((Utt_SpeechEnd > Utt_Start_l) && (ref_VAD[static_cast<int>
              (Utt_SpeechEnd) - 1] <= 0.0)) {
        Utt_SpeechEnd--;
      }

      Utt_SpeechEnd++;

      // find SpeechEnd for each utterance
      //      fprintf( 'Utt_Len is %d\n', Utt_Len);
      if (Utt_SpeechEnd - last_end >= 200.0) {
        split_align(ref_data, ref_Nsamples, ref_logVAD, deg_data, deg_Nsamples,
                    deg_logVAD, glb->Utt_Start[Utt_id], last_end, Utt_SpeechEnd,
                    glb->Utt_End[Utt_id], glb->Utt_DelayEst[Utt_id],
                    glb->Utt_DelayConf[Utt_id], glb);

        //          fprintf( '\nBest_ED1, glb.Best_D1, glb.Best_DC1 is %d, %d, %f\n',... 
        //                glb.Best_ED1, glb.Best_D1, glb.Best_DC1);
        //          fprintf( 'glb.Best_ED2, glb.Best_D2, glb.Best_DC2 is %d, %d, %f\n',... 
        //                glb.Best_ED2, glb.Best_D2, glb.Best_DC2);
        //          fprintf( 'glb.Best_BP is %d\n', glb.Best_BP);
        if ((glb->Best_DC1 > this_start) && (glb->Best_DC2 > this_start)) {
          i = static_cast<int>(((static_cast<double>(i) + 1.0) + (-1.0 -
            glb->Nutterances)) / -1.0);
          for (Utt_num = 0; Utt_num < i; Utt_num++) {
            this_start = glb->Nutterances + -static_cast<double>(Utt_num);
            speech_flag = static_cast<int>(this_start) - 1;
            count = static_cast<int>(this_start + 1.0) - 1;
            glb->Utt_DelayEst[count] = glb->Utt_DelayEst[speech_flag];
            glb->Utt_Delay[count] = glb->Utt_Delay[speech_flag];
            glb->Utt_DelayConf[count] = glb->Utt_DelayConf[speech_flag];
            glb->Utt_Start[count] = glb->Utt_Start[speech_flag];
            glb->Utt_End[count] = glb->Utt_End[speech_flag];
            glb->UttSearch_Start[count] = glb->Utt_Start[speech_flag];
            glb->UttSearch_End[count] = glb->Utt_End[speech_flag];
          }

          glb->Nutterances++;
          glb->Utt_DelayEst[Utt_id] = glb->Best_ED1;
          glb->Utt_Delay[Utt_id] = glb->Best_D1;
          glb->Utt_DelayConf[Utt_id] = glb->Best_DC1;
          glb->Utt_DelayEst[Utt_id + 1] = glb->Best_ED2;
          glb->Utt_Delay[Utt_id + 1] = glb->Best_D2;
          glb->Utt_DelayConf[Utt_id + 1] = glb->Best_DC2;
          glb->UttSearch_Start[Utt_id + 1] = glb->UttSearch_Start[Utt_id];
          glb->UttSearch_End[Utt_id + 1] = glb->UttSearch_End[Utt_id];
          if (glb->Best_D2 < glb->Best_D1) {
            glb->Utt_Start[Utt_id] = Utt_Start_l;
            glb->Utt_End[Utt_id] = glb->Best_BP;
            glb->Utt_Start[Utt_id + 1] = glb->Best_BP;
            glb->Utt_End[Utt_id + 1] = Utt_End_l;
          } else {
            glb->Utt_Start[Utt_id] = Utt_Start_l;
            glb->Utt_End[Utt_id] = glb->Best_BP + std::floor((glb->Best_D2 -
              glb->Best_D1) / (2.0 * glb->Downsample));
            glb->Utt_Start[Utt_id + 1] = glb->Best_BP - std::floor((glb->Best_D2
              - glb->Best_D1) / (2.0 * glb->Downsample));
            glb->Utt_End[Utt_id + 1] = Utt_End_l;
          }

          if ((((glb->Utt_Start[Utt_id] - 75.0) - 1.0) * glb->Downsample + 1.0)
              + glb->Best_D1 < 0.0) {
            glb->Utt_Start[Utt_id] = std::floor(((glb->Downsample - 1.0) -
              glb->Best_D1) / glb->Downsample) + 76.0;
          }

          if (((glb->Utt_End[Utt_id + 1] - 1.0) * glb->Downsample + 1.0) +
              glb->Best_D2 > deg_Nsamples - 75.0 * glb->Downsample) {
            glb->Utt_End[Utt_id + 1] = (std::floor((deg_Nsamples - glb->Best_D2)
              / glb->Downsample) - 75.0) + 1.0;
          }
        } else {
          Utt_id++;
        }
      } else {
        Utt_id++;
      }
    } else {
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  x.set_size(1, 50);
  for (i = 0; i < 50; i++) {
    this_start = glb->Utt_End[i] - glb->Utt_Start[i];
    varargin_1[i] = this_start;
    x[i] = this_start;
  }

  if (!rtIsNaN(x[0])) {
    speech_flag = 1;
  } else {
    speech_flag = 0;
    count = 2;
    exitg1 = false;
    while ((!exitg1) && (count <= 50)) {
      if (!rtIsNaN(x[count - 1])) {
        speech_flag = count;
        exitg1 = true;
      } else {
        count++;
      }
    }
  }

  if (speech_flag == 0) {
    glb->Largest_uttsize = varargin_1[0];
  } else {
    this_start = varargin_1[speech_flag - 1];
    i = speech_flag + 1;
    for (count = i; count < 51; count++) {
      Utt_Start_l = varargin_1[count - 1];
      if (this_start < Utt_Start_l) {
        this_start = Utt_Start_l;
      }
    }

    glb->Largest_uttsize = this_start;
  }

  //  fid= fopen( 'uttinfo_mat.txt', 'wt');
  //  fprintf( fid, 'Number of Utterances is:\n');
  //  fprintf( fid, '%d\n', glb.Nutterances);
  //  fprintf( fid, 'Utterance Delay Estimation:\n');
  //  fprintf( fid, '%d\n', glb.Utt_DelayEst( 1: glb.Nutterances) );
  //  fprintf( fid, 'Utterance Delay:\n');
  //  fprintf( fid, '%d\n', glb.Utt_Delay( 1: glb.Nutterances));
  //  fprintf( fid, 'Utterance Delay Confidence:\n');
  //  fprintf( fid, '%f\n', glb.Utt_DelayConf( 1: glb.Nutterances));
  //  fprintf( fid, 'Utterance Start:\n');
  //  fprintf( fid, '%d\n', glb.Utt_Start( 1: glb.Nutterances));
  //  fprintf( fid, 'Utterance End:\n');
  //  fprintf( fid, '%d\n', glb.Utt_End( 1: glb.Nutterances));
  //  fprintf( fid, 'Largest utterance length:\n');
  //  fprintf( fid, '%d\n', glb.Largest_uttsize);
  //  fclose( fid);
  //  EOF
  //  function logdata(varname,var)
  //
  //  eval(strcat(mfilename,'_',varname,'=var;'));
  //
  //  save(strcat('./test_data/',mfilename,'_',varname),strcat(mfilename,'_',varname)); 
}

void pesq_original_cpp(const coder::array<double, 1U> &ref_data, double
  ref_sampling_rate, const coder::array<double, 1U> &deg_data, double
  scores_data[], int scores_size[2])
{
  int mode_size[2];
  int i;
  char mode_data[10];
  static const char t0_f5[10] = { 'n', 'a', 'r', 'r', 'o', 'w', 'b', 'a', 'n',
    'd' };

  static const char t1_f3[8] = { 'w', 'i', 'd', 'e', 'b', 'a', 'n', 'd' };

  double glb_WB_InIIR_Hsos[5];
  static const double dv[5] = { 2.6657628, -5.3315255, 2.6657628, -1.8890331,
    0.89487434 };

  static const double dv1[5] = { 2.740826, -5.4816519, 2.740826, -1.9444777,
    0.94597794 };

  int glb_Downsample;
  int glb_InIIR_Hsos_size_idx_0;
  double glb_InIIR_Hsos_data[60];
  static const double InIIR_Hsos_16k[60] = { 0.325631521, 0.403961804,
    4.736162769, 0.365373469, 0.884811506, 0.723593055, 1.644910855, 0.633692689,
    1.032763031, 1.001616361, 0.752472096, 1.023700575, -0.08678286,
    -0.556985881, 3.287251046, 0.0, 0.0, -1.447186099, -1.817280902,
    -0.284644314, 0.268428979, -0.823749013, -0.37538899, 0.001661628,
    -0.238848661, 0.153024077, 1.753289019, 0.0, 0.0, 0.723593044, 1.249658063,
    -0.319789663, 0.602913323, 0.439731942, 0.188977609, 0.52128424, -1.07941649,
    -0.415115835, -1.859599046, -0.634626531, -0.256725271, -1.129587469,
    -1.778403899, 0.0, 0.0, -0.885778255, -0.077258216, -0.183867259,
    0.434583902, 0.696590244, 0.876284034, 0.0, 0.141536777, 0.657232737,
    0.801724355, 0.0, 0.0, 0.0, 0.247230734, 0.354324187 };

  static const double InIIR_Hsos_8k[40] = { 0.885535424, 0.895092588, 4.04952794,
    0.500002353, 0.565002834, 2.115237288, 0.912224584, 0.444617727,
    -0.885535424, 1.292907193, -7.865190042, -0.500002353, -0.241585934,
    0.919935084, -0.224397719, -0.307589321, 0.0, 0.449260174, 3.815662102, 0.0,
    -0.306009671, 1.141240051, -0.641121413, 0.141638062, -0.771070709,
    1.268869037, -1.746859852, 0.0, 0.259688659, -1.587313419, -0.246029464,
    -0.996391149, 0.0, 0.442025372, 0.786305963, 0.0, 0.249979657, 0.665935315,
    -0.55672059, 0.502251622 };

  int glb_InIIR_Nsos;
  int Align_Nfft_;
  int glb_Fs;
  int glb_Nb;
  double glb_Sp;
  int c_glb_nr_of_hz_bands_per_bark_b;
  int c_glb_centre_of_band_bark_size_;
  int c_glb_centre_of_band_hz_size_id;
  int c_glb_width_of_band_bark_size_i;
  int glb_width_of_band_hz_size_idx_1;
  int c_glb_pow_dens_correction_facto;
  int glb_abs_thresh_power_size_idx_1;
  double glb_centre_of_band_bark_data[49];
  static const double centre_of_band_bark_16k[49] = { 0.078672, 0.316341,
    0.636559, 0.961246, 1.29045, 1.624217, 1.962597, 2.305636, 2.653383,
    3.005889, 3.363201, 3.725371, 4.092449, 4.464486, 4.841533, 5.223642,
    5.610866, 6.003256, 6.400869, 6.803755, 7.211971, 7.625571, 8.044611,
    8.469146, 8.899232, 9.334927, 9.776288, 10.223374, 10.676242, 11.134952,
    11.599563, 12.070135, 12.546731, 13.029408, 13.518232, 14.013264, 14.514566,
    15.022202, 15.536238, 16.056736, 16.583761, 17.117382, 17.657663, 18.204674,
    18.758478, 19.319147, 19.886751, 20.461355, 21.043034 };

  static const double centre_of_band_bark_8k[42] = { 0.078672, 0.316341,
    0.636559, 0.961246, 1.29045, 1.624217, 1.962597, 2.305636, 2.653383,
    3.005889, 3.363201, 3.725371, 4.092449, 4.464486, 4.841533, 5.223642,
    5.610866, 6.003256, 6.400869, 6.803755, 7.211971, 7.625571, 8.044611,
    8.469146, 8.899232, 9.334927, 9.776288, 10.223374, 10.676242, 11.134952,
    11.599563, 12.070135, 12.546731, 13.029408, 13.518232, 14.013264, 14.514566,
    15.022202, 15.536238, 16.056736, 16.583761, 17.117382 };

  double glb_centre_of_band_hz_data[49];
  static const double centre_of_band_hz_16k[49] = { 7.867213, 31.634144,
    63.655895, 96.124611, 129.044968, 162.421738, 196.259659, 230.563568,
    265.338348, 300.588867, 336.320129, 372.53714, 409.244934, 446.448578,
    484.568604, 526.600586, 570.303833, 619.42334, 672.121643, 728.525696,
    785.675964, 846.835693, 909.69165, 977.063293, 1049.861694, 1129.635986,
    1217.257568, 1312.109497, 1412.501465, 1517.99939, 1628.894165, 1746.194336,
    1871.568848, 2008.776123, 2158.979248, 2326.743164, 2513.787109, 2722.48877,
    2952.58667, 3205.835449, 3492.679932, 3820.219238, 4193.938477, 4619.846191,
    5100.437012, 5636.199219, 6234.313477, 6946.734863, 7796.473633 };

  static const double centre_of_band_hz_8k[42] = { 7.867213, 31.634144,
    63.655895, 96.124611, 129.044968, 162.421738, 196.259659, 230.563568,
    265.338348, 300.588867, 336.320129, 372.53714, 409.244934, 446.448578,
    484.568604, 526.600586, 570.303833, 619.42334, 672.121643, 728.525696,
    785.675964, 846.835693, 909.69165, 977.063293, 1049.861694, 1129.635986,
    1217.257568, 1312.109497, 1412.501465, 1517.99939, 1628.894165, 1746.194336,
    1871.568848, 2008.776123, 2158.979248, 2326.743164, 2513.787109, 2722.48877,
    2952.58667, 3205.835449, 3492.679932, 3820.219238 };

  double glb_width_of_band_bark_data[49];
  static const double width_of_band_bark_16k[49] = { 0.157344, 0.317994,
    0.322441, 0.326934, 0.331474, 0.336061, 0.340697, 0.345381, 0.350114,
    0.354897, 0.359729, 0.364611, 0.369544, 0.374529, 0.379565, 0.384653,
    0.389794, 0.394989, 0.400236, 0.405538, 0.410894, 0.416306, 0.421773,
    0.427297, 0.432877, 0.438514, 0.444209, 0.449962, 0.455774, 0.461645,
    0.467577, 0.473569, 0.479621, 0.485736, 0.491912, 0.498151, 0.504454,
    0.510819, 0.51725, 0.523745, 0.530308, 0.536934, 0.543629, 0.55039, 0.55722,
    0.564119, 0.571085, 0.578125, 0.585232 };

  static const double width_of_band_bark_8k[42] = { 0.157344, 0.317994, 0.322441,
    0.326934, 0.331474, 0.336061, 0.340697, 0.345381, 0.350114, 0.354897,
    0.359729, 0.364611, 0.369544, 0.374529, 0.379565, 0.384653, 0.389794,
    0.394989, 0.400236, 0.405538, 0.410894, 0.416306, 0.421773, 0.427297,
    0.432877, 0.438514, 0.444209, 0.449962, 0.455774, 0.461645, 0.467577,
    0.473569, 0.479621, 0.485736, 0.491912, 0.498151, 0.504454, 0.510819,
    0.51725, 0.523745, 0.530308, 0.536934 };

  double glb_width_of_band_hz_data[49];
  static const double width_of_band_hz_16k[49] = { 15.734426, 31.799433,
    32.244064, 32.693359, 33.147385, 33.60614, 34.069702, 34.538116, 35.011429,
    35.489655, 35.97287, 36.461121, 36.954407, 37.452911, 40.269653, 42.311859,
    45.992554, 51.348511, 55.040527, 56.775208, 58.699402, 62.445862, 64.820923,
    69.195374, 76.745667, 84.016235, 90.825684, 97.931152, 103.348877, 107.80188,
    113.552246, 121.490601, 130.42041, 143.431763, 158.486816, 176.872803,
    198.314697, 219.549561, 240.600098, 268.702393, 306.060059, 349.937012,
    398.686279, 454.713867, 506.841797, 564.86377, 637.26123, 794.717285,
    931.068359 };

  static const double width_of_band_hz_8k[42] = { 15.734426, 31.799433,
    32.244064, 32.693359, 33.147385, 33.60614, 34.069702, 34.538116, 35.011429,
    35.489655, 35.97287, 36.461121, 36.954407, 37.452911, 40.269653, 42.311859,
    45.992554, 51.348511, 55.040527, 56.775208, 58.699402, 62.445862, 64.820923,
    69.195374, 76.745667, 84.016235, 90.825684, 97.931152, 103.348877, 107.80188,
    113.552246, 121.490601, 130.42041, 143.431763, 158.486816, 176.872803,
    198.314697, 219.549561, 240.600098, 268.702393, 306.060059, 349.937012 };

  double d_glb_pow_dens_correction_facto[49];
  static const double pow_dens_correction_factor_16k[49] = { 100.0, 99.999992,
    100.0, 100.000008, 100.000008, 100.000015, 99.999992, 99.999969, 50.000027,
    100.0, 99.999969, 100.000015, 99.999947, 100.000061, 53.047077, 110.000046,
    117.991989, 65.0, 68.760147, 69.999931, 71.428818, 75.000038, 76.843384,
    80.968781, 88.646126, 63.864388, 68.15535, 72.547775, 75.584831, 58.379192,
    80.950836, 64.135651, 54.384785, 73.821884, 64.437073, 59.176456, 65.521278,
    61.399822, 58.144047, 57.004543, 64.126297, 54.311001, 61.114979, 55.077751,
    56.849335, 55.628868, 53.137054, 54.985844, 79.546974 };

  static const double pow_dens_correction_factor_8k[42] = { 100.0, 99.999992,
    100.0, 100.000008, 100.000008, 100.000015, 99.999992, 99.999969, 50.000027,
    100.0, 99.999969, 100.000015, 99.999947, 100.000061, 53.047077, 110.000046,
    117.991989, 65.0, 68.760147, 69.999931, 71.428818, 75.000038, 76.843384,
    80.968781, 88.646126, 63.864388, 68.15535, 72.547775, 75.584831, 58.379192,
    80.950836, 64.135651, 54.384785, 73.821884, 64.437073, 59.176456, 65.521278,
    61.399822, 58.144047, 57.004543, 64.126297, 59.248363 };

  double glb_abs_thresh_power_data[49];
  static const double abs_thresh_power_16k[49] = { 5.1286152E+7, 2.4547095E+6,
    70794.59375, 4897.788574, 1174.897705, 389.045166, 104.71286, 45.70882,
    17.782795, 9.772372, 4.897789, 3.090296, 1.905461, 1.258925, 0.977237,
    0.724436, 0.562341, 0.457088, 0.389045, 0.331131, 0.295121, 0.269153,
    0.25704, 0.251189, 0.251189, 0.251189, 0.251189, 0.263027, 0.288403, 0.30903,
    0.338844, 0.371535, 0.398107, 0.436516, 0.467735, 0.489779, 0.501187,
    0.501187, 0.512861, 0.524807, 0.524807, 0.524807, 0.512861, 0.47863, 0.42658,
    0.371535, 0.363078, 0.416869, 0.537032 };

  static const double abs_thresh_power_8k[42] = { 5.1286152E+7, 2.4547095E+6,
    70794.59375, 4897.788574, 1174.897705, 389.045166, 104.71286, 45.70882,
    17.782795, 9.772372, 4.897789, 3.090296, 1.905461, 1.258925, 0.977237,
    0.724436, 0.562341, 0.457088, 0.389045, 0.331131, 0.295121, 0.269153,
    0.25704, 0.251189, 0.251189, 0.251189, 0.251189, 0.263027, 0.288403, 0.30903,
    0.338844, 0.371535, 0.398107, 0.436516, 0.467735, 0.489779, 0.501187,
    0.501187, 0.512861, 0.524807, 0.524807, 0.524807 };

  signed char d_glb_nr_of_hz_bands_per_bark_b[49];
  static const signed char c_nr_of_hz_bands_per_bark_band_[49] = { 1, 1, 1, 1, 1,
    1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 3,
    4, 5, 4, 5, 6, 6, 7, 8, 9, 9, 12, 12, 15, 16, 18, 21, 25, 20 };

  static const signed char nr_of_hz_bands_per_bark_band_8k[42] = { 1, 1, 1, 1, 1,
    1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 3,
    4, 5, 4, 5, 6, 6, 7, 8, 9, 9, 11 };

  int count;
  coder::array<double, 2U> b_ref_data;
  double glb_Window_data[1024];
  unsigned int ref_Nsamples;
  int unnamed_idx_1;
  int unnamed_idx_1_tmp;
  coder::array<double, 2U> c_ref_data;
  coder::array<double, 2U> b_deg_data;
  unsigned int deg_Nsamples;
  unsigned int maxNsamples;
  b_struct_T expl_temp;
  struct_T glb;
  coder::array<double, 2U> c_deg_data;
  boolean_T b_bool;
  char switch_expression_tmp_data[10];
  static const char cv[128] = { '\x00', '\x01', '\x02', '\x03', '\x04', '\x05',
    '\x06', '\x07', '\x08', '\x09', '\x0a', '\x0b', '\x0c', '\x0d', '\x0e',
    '\x0f', '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17',
    '\x18', '\x19', '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', ' ', '!',
    '\"', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0',
    '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
    '@', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
    'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '[', '\\', ']',
    '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{',
    '|', '}', '~', '\x7f' };

  int exitg1;
  coder::array<char, 2U> b_mode_data;
  static const char t1_f4[9] = { '+', 'w', 'i', 'd', 'e', 'b', 'a', 'n', 'd' };

  coder::array<double, 2U> ref_VAD;
  coder::array<double, 2U> ref_logVAD;
  coder::array<double, 2U> deg_VAD;
  coder::array<double, 2U> deg_logVAD;
  c_struct_T b_glb;
  coder::array<char, 2U> c_mode_data;

  //  ----------------------------------------------------------------------
  //             PESQ objective speech quality measure
  //             (narrowband and wideband implementations)
  //
  //    This function implements the PESQ measure based on the ITU standards
  //    P.862 [1] and P.862.1 [2] for narrowband speech and P.862.2 for
  //    wideband speech [3].
  //
  //
  //    Usage:  scores = pesq( cleanFile, enhancedFile )
  //
  //          cleanFile     - clean input file in .wav format sampled at
  //                          sampling frequency Fs=8 kHz or Fs=16 kHz
  //                          for narrowband or wideband assessment,
  //                          respectively.
  //
  //          enhancedFile  - enhanced output file in .wav format sampled
  //                          at same sampling frequency as the cleanFile
  //
  //          scores        - For narrowband speech, two scores are returned,
  //                          one for the raw PESQ value [1] (first value) and
  //                          one for the MOS-mapped score value [2] (second value). 
  //                          For wideband speech, only the MOS-mapped value
  //                          is returned [3].
  //
  //   Example call:  scores = pesq('sp04.wav', 'enhanced.wav')
  //
  //
  //   References:
  //
  //    [1] ITU (2000). Perceptual evaluation of speech quality (PESQ), and
  //        objective method for end-to-end speech quality assessment of
  //        narrowband telephone networks and speech codecs. ITU-T
  //        Recommendation P.862
  //
  //    [2] ITU (2003).  Mapping function for transforming P.862 raw result
  //        scores to MOS-LQO, ITU-T Recommendation P. 862.1
  //
  //    [3] ITU (2007). Wideband extension to Recommendation P.862 for the
  //        assessment of wideband telephone networks and speech codecs. ITU-T
  //        Recommendation P.862.2
  //
  //
  //    Authors: Yi Hu, Kamil Wojcicki and Philipos C. Loizou
  //
  //
  //  Copyright (c) 2006, 2012 by Philipos C. Loizou
  //  $Revision: 2.0 $  $Date: 5/14/2012 $
  //  ----------------------------------------------------------------------
  //      if nargin ==0, fprintf('Usage: pesq( ref_wav, deg_wav )\n');
  //                     fprintf('       ref_wav = reference input filename\n'); 
  //                     fprintf('       deg_wav = degraded output filename\n\n'); 
  //                     fprintf('For more help, type: help pesq\n\n');
  //                     return;
  //      elseif nargin > 2, error('%s.m: Incorrect number of input arguments.\nFor usage help type: help %s',mfilename,mfilename); 
  //      end
  //      if ~isstr( ref_wav ), error( '%s.m: First input argumnet has to be a reference wav filename as string.\nFor usage help type: help %s',mfilename,mfilename); end; 
  //      if ~isstr( deg_wav ), error( '%s.m: Second input argumnet has to be a processed wav filename as string.\nFor usage help type: help %s',mfilename,mfilename); end; 
  //
  //      if ~exist( ref_wav, 'file' ), error( '%s.m: Reference wav file: %s not found.',mfilename,ref_wav); end; 
  //      if ~exist( deg_wav, 'file' ), error( '%s.m: Processed wav file: %s not found.',mfilename,deg_wav); end; 
  //      [ ref_data, ref_sampling_rate ] = audioread( ref_wav );
  //      [ deg_data, deg_sampling_rate ] = audioread( deg_wav );
  if (ref_sampling_rate == 8000.0) {
    mode_size[0] = 1;
    mode_size[1] = 10;
    for (i = 0; i < 10; i++) {
      mode_data[i] = t0_f5[i];
    }
  } else {
    if (ref_sampling_rate == 16000.0) {
      mode_size[0] = 1;
      mode_size[1] = 8;
      for (i = 0; i < 8; i++) {
        mode_data[i] = t1_f3[i];
      }
    }
  }

  //      clearvars -global Downsample DATAPADDING_MSECS SEARCHBUFFER Fs WHOLE_SIGNAL Align_Nfft Window 
  //      global Downsample DATAPADDING_MSECS SEARCHBUFFER Fs WHOLE_SIGNAL
  //      global Align_Nfft Window
  //      global glb.Downsample glb.InIIR_Hsos glb.InIIR_Nsos glb.Align_Nfft
  //      global glb.DATAPADDING_MSEC glb.SEARCHBUFFER glb.Fs glb.glb.MINSPEECHLGTH glb.JOINSPEECHLGTH 
  //
  //      global glb.Nutterances glb.Largest_uttsize glb.Nsurf_samples glb.Crude_DelayEst 
  //      global glb.Crude_DelayConf glb.UttSearch_Start glb.UttSearch_End glb.Utt_DelayEst 
  //      global glb.Utt_Delay glb.Utt_DelayConf glb.Utt_Start glb.Utt_End
  //      global glb.MAXNUTTERANCES glb.WHOLE_SIGNAL
  //      global glb.pesq_mos glb.subj_mos glb.cond_nr glb.MINUTTLENGTH
  //      global glb.CALIBRATE glb.Nfmax glb.Nb glb.Sl glb.Sp
  //      global glb.nr_of_hz_bands_per_bark_band glb.centre_of_band_bark
  //      global glb.width_of_band_hz glb.centre_of_band_hz glb.width_of_band_bark  
  //      global glb.pow_dens_correction_factor glb.abs_thresh_power
  //      glb = struct();
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  // zeros( 1, glb.MAXNUTTERANCES);
  //  KKW ---------
  //      global glb.WB_InIIR_Nsos glb.WB_InIIR_Hsos
  switch (static_cast<int>(ref_sampling_rate)) {
   case 8000:
    for (i = 0; i < 5; i++) {
      glb_WB_InIIR_Hsos[i] = dv[i];
    }
    break;

   case 16000:
    for (i = 0; i < 5; i++) {
      glb_WB_InIIR_Hsos[i] = dv1[i];
    }
    break;
  }

  //  -------------
  if (ref_sampling_rate == 16000.0) {
    glb_Downsample = 64;
    glb_InIIR_Hsos_size_idx_0 = 12;
    std::memcpy(&glb_InIIR_Hsos_data[0], &InIIR_Hsos_16k[0], 60U * sizeof(double));
    glb_InIIR_Nsos = 12;
    Align_Nfft_ = 1024;
    glb_Fs = 16000;
    glb_Nb = 49;
    glb_Sp = 6.910853E-6;
    c_glb_nr_of_hz_bands_per_bark_b = 49;
    c_glb_centre_of_band_bark_size_ = 49;
    c_glb_centre_of_band_hz_size_id = 49;
    c_glb_width_of_band_bark_size_i = 49;
    glb_width_of_band_hz_size_idx_1 = 49;
    c_glb_pow_dens_correction_facto = 49;
    glb_abs_thresh_power_size_idx_1 = 49;
    std::memcpy(&glb_centre_of_band_bark_data[0], &centre_of_band_bark_16k[0],
                49U * sizeof(double));
    std::memcpy(&glb_centre_of_band_hz_data[0], &centre_of_band_hz_16k[0], 49U *
                sizeof(double));
    std::memcpy(&glb_width_of_band_bark_data[0], &width_of_band_bark_16k[0], 49U
                * sizeof(double));
    std::memcpy(&glb_width_of_band_hz_data[0], &width_of_band_hz_16k[0], 49U *
                sizeof(double));
    std::memcpy(&d_glb_pow_dens_correction_facto[0],
                &pow_dens_correction_factor_16k[0], 49U * sizeof(double));
    std::memcpy(&glb_abs_thresh_power_data[0], &abs_thresh_power_16k[0], 49U *
                sizeof(double));
    for (i = 0; i < 49; i++) {
      d_glb_nr_of_hz_bands_per_bark_b[i] = c_nr_of_hz_bands_per_bark_band_[i];
    }
  } else {
    // (sampling_rate== fs_8k)
    glb_Downsample = 32;
    glb_InIIR_Hsos_size_idx_0 = 8;
    std::memcpy(&glb_InIIR_Hsos_data[0], &InIIR_Hsos_8k[0], 40U * sizeof(double));
    glb_InIIR_Nsos = 8;
    Align_Nfft_ = 512;
    glb_Fs = 8000;
    glb_Nb = 42;
    glb_Sp = 2.764344E-5;
    c_glb_nr_of_hz_bands_per_bark_b = 42;
    c_glb_centre_of_band_bark_size_ = 42;
    c_glb_centre_of_band_hz_size_id = 42;
    c_glb_width_of_band_bark_size_i = 42;
    glb_width_of_band_hz_size_idx_1 = 42;
    c_glb_pow_dens_correction_facto = 42;
    glb_abs_thresh_power_size_idx_1 = 42;
    std::memcpy(&glb_centre_of_band_bark_data[0], &centre_of_band_bark_8k[0],
                42U * sizeof(double));
    std::memcpy(&glb_centre_of_band_hz_data[0], &centre_of_band_hz_8k[0], 42U *
                sizeof(double));
    std::memcpy(&glb_width_of_band_bark_data[0], &width_of_band_bark_8k[0], 42U *
                sizeof(double));
    std::memcpy(&glb_width_of_band_hz_data[0], &width_of_band_hz_8k[0], 42U *
                sizeof(double));
    std::memcpy(&d_glb_pow_dens_correction_facto[0],
                &pow_dens_correction_factor_8k[0], 42U * sizeof(double));
    std::memcpy(&glb_abs_thresh_power_data[0], &abs_thresh_power_8k[0], 42U *
                sizeof(double));
    for (i = 0; i < 42; i++) {
      d_glb_nr_of_hz_bands_per_bark_b[i] = nr_of_hz_bands_per_bark_band_8k[i];
    }
  }

  for (count = 0; count < Align_Nfft_; count++) {
    glb_Window_data[count] = 0.5 * (1.0 - std::cos(6.28318530717959 *
      static_cast<double>(count) / static_cast<double>(Align_Nfft_)));
  }

  b_ref_data.set_size(1, ref_data.size(0));
  count = ref_data.size(0);
  for (i = 0; i < count; i++) {
    b_ref_data[i] = ref_data[i] * 32768.0;
  }

  ref_Nsamples = static_cast<unsigned int>(b_ref_data.size(1)) + 150 *
    glb_Downsample;
  unnamed_idx_1 = 75 * glb_Downsample;
  unnamed_idx_1_tmp = static_cast<int>(320.0 * (static_cast<double>(glb_Fs) /
    1000.0) + 75.0 * static_cast<double>(glb_Downsample));
  c_ref_data.set_size(1, ((unnamed_idx_1 + b_ref_data.size(1)) +
    unnamed_idx_1_tmp));
  for (i = 0; i < unnamed_idx_1; i++) {
    c_ref_data[i] = 0.0;
  }

  count = b_ref_data.size(1);
  for (i = 0; i < count; i++) {
    c_ref_data[i + unnamed_idx_1] = b_ref_data[i];
  }

  for (i = 0; i < unnamed_idx_1_tmp; i++) {
    c_ref_data[(i + unnamed_idx_1) + b_ref_data.size(1)] = 0.0;
  }

  b_ref_data.set_size(1, c_ref_data.size(1));
  count = c_ref_data.size(0) * c_ref_data.size(1);
  for (i = 0; i < count; i++) {
    b_ref_data[i] = c_ref_data[i];
  }

  b_deg_data.set_size(1, deg_data.size(0));
  count = deg_data.size(0);
  for (i = 0; i < count; i++) {
    b_deg_data[i] = deg_data[i] * 32768.0;
  }

  deg_Nsamples = static_cast<unsigned int>(b_deg_data.size(1)) + 150 *
    glb_Downsample;
  unnamed_idx_1 = 75 * glb_Downsample;
  c_ref_data.set_size(1, ((unnamed_idx_1 + b_deg_data.size(1)) +
    unnamed_idx_1_tmp));
  for (i = 0; i < unnamed_idx_1; i++) {
    c_ref_data[i] = 0.0;
  }

  count = b_deg_data.size(1);
  for (i = 0; i < count; i++) {
    c_ref_data[i + unnamed_idx_1] = b_deg_data[i];
  }

  for (i = 0; i < unnamed_idx_1_tmp; i++) {
    c_ref_data[(i + unnamed_idx_1) + b_deg_data.size(1)] = 0.0;
  }

  b_deg_data.set_size(1, c_ref_data.size(1));
  count = c_ref_data.size(0) * c_ref_data.size(1);
  for (i = 0; i < count; i++) {
    b_deg_data[i] = c_ref_data[i];
  }

  if (ref_Nsamples > deg_Nsamples) {
    maxNsamples = ref_Nsamples;
  } else {
    maxNsamples = deg_Nsamples;
  }

  expl_temp.Window.size[0] = 1;
  expl_temp.Window.size[1] = Align_Nfft_;
  std::memcpy(&expl_temp.Window.data[0], &glb_Window_data[0], Align_Nfft_ *
              sizeof(double));
  expl_temp.abs_thresh_power.size[0] = 1;
  expl_temp.abs_thresh_power.size[1] = glb_abs_thresh_power_size_idx_1;
  std::memcpy(&expl_temp.abs_thresh_power.data[0], &glb_abs_thresh_power_data[0],
              glb_abs_thresh_power_size_idx_1 * sizeof(double));
  expl_temp.pow_dens_correction_factor.size[0] = 1;
  expl_temp.pow_dens_correction_factor.size[1] = c_glb_pow_dens_correction_facto;
  std::memcpy(&expl_temp.pow_dens_correction_factor.data[0],
              &d_glb_pow_dens_correction_facto[0],
              c_glb_pow_dens_correction_facto * sizeof(double));
  expl_temp.width_of_band_hz.size[0] = 1;
  expl_temp.width_of_band_hz.size[1] = glb_width_of_band_hz_size_idx_1;
  std::memcpy(&expl_temp.width_of_band_hz.data[0], &glb_width_of_band_hz_data[0],
              glb_width_of_band_hz_size_idx_1 * sizeof(double));
  expl_temp.width_of_band_bark.size[0] = 1;
  expl_temp.width_of_band_bark.size[1] = c_glb_width_of_band_bark_size_i;
  std::memcpy(&expl_temp.width_of_band_bark.data[0],
              &glb_width_of_band_bark_data[0], c_glb_width_of_band_bark_size_i *
              sizeof(double));
  expl_temp.centre_of_band_hz.size[0] = 1;
  expl_temp.centre_of_band_hz.size[1] = c_glb_centre_of_band_hz_size_id;
  std::memcpy(&expl_temp.centre_of_band_hz.data[0], &glb_centre_of_band_hz_data
              [0], c_glb_centre_of_band_hz_size_id * sizeof(double));
  expl_temp.centre_of_band_bark.size[0] = 1;
  expl_temp.centre_of_band_bark.size[1] = c_glb_centre_of_band_bark_size_;
  std::memcpy(&expl_temp.centre_of_band_bark.data[0],
              &glb_centre_of_band_bark_data[0], c_glb_centre_of_band_bark_size_ *
              sizeof(double));
  expl_temp.nr_of_hz_bands_per_bark_band.size[0] = 1;
  expl_temp.nr_of_hz_bands_per_bark_band.size[1] =
    c_glb_nr_of_hz_bands_per_bark_b;
  for (i = 0; i < c_glb_nr_of_hz_bands_per_bark_b; i++) {
    expl_temp.nr_of_hz_bands_per_bark_band.data[i] =
      d_glb_nr_of_hz_bands_per_bark_b[i];
  }

  expl_temp.Align_Nfft = Align_Nfft_;
  expl_temp.InIIR_Nsos = glb_InIIR_Nsos;
  expl_temp.InIIR_Hsos.size[0] = glb_InIIR_Hsos_size_idx_0;
  expl_temp.InIIR_Hsos.size[1] = 5;
  count = glb_InIIR_Hsos_size_idx_0 * 5;
  std::memcpy(&expl_temp.InIIR_Hsos.data[0], &glb_InIIR_Hsos_data[0], count *
              sizeof(double));
  expl_temp.Fs = glb_Fs;
  expl_temp.Sp = glb_Sp;
  expl_temp.Sl = 0.1866055;
  expl_temp.Nb = glb_Nb;
  expl_temp.Downsample = glb_Downsample;
  for (i = 0; i < 5; i++) {
    expl_temp.WB_InIIR_Hsos[i] = glb_WB_InIIR_Hsos[i];
  }

  expl_temp.WB_InIIR_Nsos = 1.0;
  expl_temp.Best_BP = 0.0;
  expl_temp.Best_D2 = 0.0;
  expl_temp.Best_ED2 = 0.0;
  expl_temp.Best_D1 = 0.0;
  expl_temp.Best_ED1 = 0.0;
  expl_temp.Best_DC2 = 0.0;
  expl_temp.Best_DC1 = 0.0;
  expl_temp.Largest_uttsize = 0.0;
  expl_temp.Nutterances = 1.0;
  expl_temp.Crude_DelayConf = 0.0;
  expl_temp.Crude_DelayEst = 0.0;
  expl_temp.JOINSPEECHLGTH = 50.0;
  expl_temp.MINSPEECHLGTH = 4.0;
  expl_temp.SEARCHBUFFER = 75.0;
  expl_temp.DATAPADDING_MSEC = 320.0;
  std::memset(&expl_temp.Utt_End[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.Utt_Start[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.Utt_DelayConf[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.Utt_Delay[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.Utt_DelayEst[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.UttSearch_End[0], 0, 50U * sizeof(double));
  std::memset(&expl_temp.UttSearch_Start[0], 0, 50U * sizeof(double));
  expl_temp.WHOLE_SIGNAL = -1.0;
  expl_temp.MINUTTLENGTH = 50.0;
  expl_temp.MAXNUTTERANCES = 50.0;
  expl_temp.Nfmax = 512.0;
  expl_temp.CALIBRATE = 0.0;
  fix_power_level(b_ref_data, static_cast<double>(ref_Nsamples), static_cast<
                  double>(maxNsamples), &expl_temp, c_ref_data, &glb);
  b_ref_data.set_size(1, c_ref_data.size(1));
  count = c_ref_data.size(0) * c_ref_data.size(1);
  for (i = 0; i < count; i++) {
    b_ref_data[i] = c_ref_data[i];
  }

  b_fix_power_level(b_deg_data, static_cast<double>(deg_Nsamples), static_cast<
                    double>(maxNsamples), &glb, c_deg_data);
  b_deg_data.set_size(1, c_deg_data.size(1));
  count = c_deg_data.size(0) * c_deg_data.size(1);
  for (i = 0; i < count; i++) {
    b_deg_data[i] = c_deg_data[i];
  }

  //  KKW ---------
  i = mode_size[1];
  for (count = 0; count < i; count++) {
    switch_expression_tmp_data[count] = cv[static_cast<int>(mode_data[count])];
  }

  b_bool = false;
  if (mode_size[1] == 10) {
    count = 0;
    do {
      exitg1 = 0;
      if (count < 10) {
        if (switch_expression_tmp_data[count] != t0_f5[count]) {
          exitg1 = 1;
        } else {
          count++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    count = 0;
  } else {
    b_bool = false;
    if (mode_size[1] == 8) {
      count = 0;
      do {
        exitg1 = 0;
        if (count < 8) {
          if (switch_expression_tmp_data[count] != t1_f3[count]) {
            exitg1 = 1;
          } else {
            count++;
          }
        } else {
          b_bool = true;
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    if (b_bool) {
      count = 1;
    } else {
      b_bool = false;
      if (mode_size[1] == 9) {
        count = 0;
        do {
          exitg1 = 0;
          if (count < 9) {
            if (switch_expression_tmp_data[count] != t1_f4[count]) {
              exitg1 = 1;
            } else {
              count++;
            }
          } else {
            b_bool = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }

      if (b_bool) {
        count = 1;
      } else {
        count = -1;
      }
    }
  }

  switch (count) {
   case 0:
    b_apply_filter(c_ref_data, static_cast<double>(ref_Nsamples), &glb,
                   b_ref_data);
    b_apply_filter(c_deg_data, static_cast<double>(deg_Nsamples), &glb,
                   b_deg_data);

    //              logdata('apply_filter_ref_data',ref_data);
    //              logdata('apply_filter_deg_data',deg_data);
    break;

   case 1:
    apply_filters_WB(c_ref_data, &glb, b_ref_data);
    apply_filters_WB(c_deg_data, &glb, b_deg_data);
    break;

   default:
    b_sprintf(mode_data, mode_size, b_mode_data);
    break;
  }

  //  -------------
  //
  //  fid= fopen( 'log_mat_ref.txt', 'wt');
  //  fprintf( fid, '%f\n', ref_data);
  //  fclose( fid);
  //
  //  fid= fopen( 'log_mat_deg.txt', 'wt');
  //  fprintf( fid, '%f\n', deg_data);
  //  fclose( fid);
  //  % to save time, read from data file ========
  //  fid= fopen( 'log_mat_ref.txt', 'rt');
  //  ref_data= fscanf( fid, '%f\n');
  //  ref_data= ref_data';
  //  fclose( fid);
  //  ref_Nsamples= length( ref_data)- glb.DATAPADDING_MSEC* (glb.Fs/ 1000);
  //
  //  fid= fopen( 'log_mat_deg.txt', 'rt');
  //  deg_data= fscanf( fid, '%f\n');
  //  deg_data= deg_data';
  //  fclose( fid);
  //  deg_Nsamples= length( deg_data)- glb.DATAPADDING_MSEC* (glb.Fs/ 1000);
  //  % the above part will be commented after debugging ========
  //  for later use in psychoacoustical model
  input_filter(b_ref_data, static_cast<double>(ref_Nsamples), b_deg_data,
               static_cast<double>(deg_Nsamples), &glb, c_ref_data, c_deg_data);

  //      logdata('ref_data_if', ref_data);
  //  fid= fopen( 'log_mat_ref_tovad.txt', 'wt');
  //  fprintf( fid, '%f\n', ref_data);
  //  fclose( fid);
  //
  //  fid= fopen( 'log_mat_deg_tovad.txt', 'wt');
  //  fprintf( fid, '%f\n', deg_data);
  //  fclose( fid);
  apply_VAD(c_ref_data, static_cast<double>(ref_Nsamples), &glb, ref_VAD,
            ref_logVAD);
  apply_VAD(c_deg_data, static_cast<double>(deg_Nsamples), &glb, deg_VAD,
            deg_logVAD);

  //  subplot( 2, 2, 1); plot( ref_VAD); title( 'ref\_VAD');
  //  subplot( 2, 2, 2); plot( ref_logVAD); title( 'ref\_logVAD');
  //
  //  subplot( 2, 2, 3); plot( deg_VAD); title( 'deg\_VAD');
  //  subplot( 2, 2, 4); plot( deg_logVAD); title( 'deg\_logVAD');
  //
  //  fid= fopen( 'mat_ref_vad.txt', 'wt');
  //  fprintf( fid, '%f\n', ref_VAD);
  //  fclose( fid);
  //
  //  fid= fopen( 'mat_ref_logvad.txt', 'wt');
  //  fprintf( fid, '%f\n', ref_logVAD);
  //  fclose( fid);
  //
  //  fid= fopen( 'mat_deg_vad.txt', 'wt');
  //  fprintf( fid, '%f\n', deg_VAD);
  //  fclose( fid);
  //
  //  fid= fopen( 'mat_deg_logvad.txt', 'wt');
  //  fprintf( fid, '%f\n', deg_logVAD);
  //  fclose( fid);
  //
  crude_align(ref_logVAD, static_cast<double>(ref_Nsamples), deg_logVAD,
              static_cast<double>(deg_Nsamples), -1.0, &glb);
  utterance_locate(c_ref_data, static_cast<double>(ref_Nsamples), ref_VAD,
                   ref_logVAD, c_deg_data, static_cast<double>(deg_Nsamples),
                   deg_logVAD, &glb);

  //  make ref_data and deg_data equal length
  if (ref_Nsamples < deg_Nsamples) {
    b_ref_data[static_cast<int>(static_cast<double>(deg_Nsamples) + 320.0 *
      (glb.Fs / 1000.0)) - 1] = 0.0;
  } else {
    if (ref_Nsamples > deg_Nsamples) {
      b_deg_data[static_cast<int>(static_cast<double>(ref_Nsamples) + 320.0 *
        (glb.Fs / 1000.0)) - 1] = 0.0;
    }
  }

  //      logdata('ref_data_pm',ref_data);
  //      logdata('deg_data_pm',deg_data);
  pesq_psychoacoustic_model(b_ref_data, static_cast<double>(ref_Nsamples),
    b_deg_data, static_cast<double>(deg_Nsamples), &glb, &glb_Sp, &b_glb);

  //  KKW ---------
  b_bool = false;
  if (mode_size[1] == 10) {
    count = 0;
    do {
      exitg1 = 0;
      if (count < 10) {
        if (switch_expression_tmp_data[count] != t0_f5[count]) {
          exitg1 = 1;
        } else {
          count++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    count = 0;
  } else {
    b_bool = false;
    if (mode_size[1] == 8) {
      count = 0;
      do {
        exitg1 = 0;
        if (count < 8) {
          if (switch_expression_tmp_data[count] != t1_f3[count]) {
            exitg1 = 1;
          } else {
            count++;
          }
        } else {
          b_bool = true;
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    if (b_bool) {
      count = 1;
    } else {
      b_bool = false;
      if (mode_size[1] == 9) {
        count = 0;
        do {
          exitg1 = 0;
          if (count < 9) {
            if (switch_expression_tmp_data[count] != t1_f4[count]) {
              exitg1 = 1;
            } else {
              count++;
            }
          } else {
            b_bool = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }

      if (b_bool) {
        count = 1;
      } else {
        count = -1;
      }
    }
  }

  switch (count) {
   case 0:
    //  NB: P.862.1->P.800.1 (PESQ_MOS->MOS_LQO)
    scores_size[0] = 1;
    scores_size[1] = 2;
    scores_data[0] = glb_Sp;
    scores_data[1] = 3.9999999999999996 / (std::exp(-1.4945 * glb_Sp + 4.6607) +
      1.0) + 0.999;
    break;

   case 1:
    //  WB: P.862.2->P.800.1 (PESQ_MOS->MOS_LQO)
    scores_size[0] = 1;
    scores_size[1] = 1;
    scores_data[0] = 3.9999999999999996 / (std::exp(-1.3669 * glb_Sp + 3.8224) +
      1.0) + 0.999;
    break;

   default:
    scores_size[0] = 1;
    scores_size[1] = 1;
    scores_data[0] = -1.0;
    b_sprintf(mode_data, mode_size, c_mode_data);
    break;
  }
}

// End of code generation (pesq_original_cpp.cpp)
