
#include <FFTImplementationCallback.h>
#include <composite_cpp.h>
#include <composite_cpp_terminate.h>
#include <gencoswin.h>
#include <rt_nonfinite.h>
#include <composite.h>

/*!
 * \brief This overloaded () operator computes the objective ITU.835 measurements.
 *
 * \details Based on clear speech signal and an enhanced speech signal, asseses the quality
 * of the enhanced signal in three categories:
 *
 * \p sig - signal quality
 * \p bak - the background noise quality
 * \p ovl - the overal quality of the enhanced sample
 *
 * \param[in] clean_signal pointer to memory space containing samples of the clean signal
 * \param[in] clean_signal_length number of samples in the \p clean_signal
 * \param[in] sampling_rate_clean_signal the sampling rate of the \p clean_signal
 * \param[in] enhanced_signal pointer to memory space containing samples of the enhanced signal
 * \param[in] enhanced_signal_length number of samples in the \p enhanced_signal
 * \param[in] sampling_rate_enhanced_signal the sampling rate of the \p enhanced_signal
 *
 * \return structure with fields
 * \stuct
 * \var double csig signal quality
 * \var double cbak background quality
 * \var double covl overal quality
 *
 */
t2r::OutputStructure t2r::Composite::operator()(const float *clean_signal, const int clean_signal_length, const double sampling_rate_clean_signal, \
      const float *enhanced_signal, const int enhanced_signal_length, const double sampling_rate_enhanced_signal){

    t2r::OutputStructure output;

    int Csig_size[2];
    int Cbak_size[2];
    int Covl_size[2];

    coder::array<double, 1U> clean_signal_array;
    coder::array<double, 1U> enhanced_signal_array;

    clean_signal_array.set_size(clean_signal_length);

    enhanced_signal_array.set_size(enhanced_signal_length);

    for(size_t i = 0; i < clean_signal_array.size(0); i++){
      clean_signal_array[i] = (double)clean_signal[i];
    }
  
    for(size_t i = 0; i < enhanced_signal_array.size(0); i++){
      enhanced_signal_array[i] = (double)enhanced_signal[i];
    }


  // Call the entry-point 'composite_cpp'.
  composite_cpp(clean_signal_array, enhanced_signal_array, sampling_rate_clean_signal, sampling_rate_enhanced_signal, &output.csig,
                Csig_size, &output.cbak, Cbak_size, &output.covl, Covl_size);

   
  return output;
} 

