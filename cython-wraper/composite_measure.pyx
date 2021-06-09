import numpy as np
cimport numpy as np


cdef extern from "composite.h" namespace "t2r":

    cdef cppclass OutputStructure:
        double csig;
        double cbak;
        double covl;

    cdef cppclass Composite:
        Composite()
        OutputStructure operator()(const float *clean_signal, const int clean_signal_length, const double sampling_rate_clean_signal, const float *enhanced_signal, const int enhanced_signal_length, const double sampling_rate_enhanced_signal)


def composite(np.ndarray[np.float32_t, ndim=1] clean_signal, double sampling_rate_clean_signal, np.ndarray[np.float32_t, ndim=1] enhanced_signal, double sampling_rate_enhanced_signal  ):

    clean_signal = np.ascontiguousarray(clean_signal)
    enhanced_signal = np.ascontiguousarray(enhanced_signal)

    cdef OutputStructure out
    cdef Composite comp

    out = comp(&clean_signal[0], clean_signal.shape[0], sampling_rate_clean_signal, &enhanced_signal[0], enhanced_signal.shape[0], sampling_rate_enhanced_signal)
    
    return out.csig, out.cbak, out.covl


