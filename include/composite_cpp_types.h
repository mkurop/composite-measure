//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  composite_cpp_types.h
//
//  Code generation for function 'composite_cpp_types'
//


#ifndef COMPOSITE_CPP_TYPES_H
#define COMPOSITE_CPP_TYPES_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#ifdef _MSC_VER

#pragma warning(push)
#pragma warning(disable : 4251)

#endif

// Type Definitions
class FFTImplementationCallback
{
 public:
  static void get_algo_sizes(int nfft, boolean_T useRadix2, int *n2blue, int
    *nRows);
  static void dobluesteinfft(const coder::array<double, 1U> &x, int n2blue, int
    nfft, const coder::array<double, 2U> &costab, const coder::array<double, 2U>
    &sintab, const coder::array<double, 2U> &sintabinv, coder::array<creal_T, 1U>
    &y);
  static void r2br_r2dit_trig_impl(const coder::array<creal_T, 1U> &x, int
    unsigned_nRows, const coder::array<double, 2U> &costab, const coder::array<
    double, 2U> &sintab, coder::array<creal_T, 1U> &y);
  static void doHalfLengthRadix2(const coder::array<double, 1U> &x, coder::array<
    creal_T, 1U> &y, int unsigned_nRows, const coder::array<double, 2U> &costab,
    const coder::array<double, 2U> &sintab);
 protected:
  static void doHalfLengthBluestein(const coder::array<double, 1U> &x, coder::
    array<creal_T, 1U> &y, int nrowsx, int nRows, int nfft, const coder::array<
    creal_T, 1U> &wwc, const coder::array<double, 2U> &costab, const coder::
    array<double, 2U> &sintab, const coder::array<double, 2U> &costabinv, const
    coder::array<double, 2U> &sintabinv);
};

#define MAX_THREADS                    omp_get_max_threads()
#ifdef _MSC_VER

#pragma warning(pop)

#endif
#endif

// End of code generation (composite_cpp_types.h)
