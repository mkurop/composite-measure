#ifndef _API_COMPOSITE_H
#define _API_COMPOSITE_H

namespace t2r{

  struct OutputStructure{ 
    double csig{0};
    double cbak{0};
    double covl{0};
  };

class Composite{

  public:

    Composite() = default;

    OutputStructure operator()(const float *clean_signal, const int clean_signal_length, const double sampling_rate_clean_signal, const float *enhanced_signal, const int enhanced_signal_length, const double sampling_rate_enhanced_signal);

};

}

#endif
