

//***********************************************************************
// This automatically generated example C++ main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************

// Include files
#include <composite.h>
#include "main.h"
#include "FFTImplementationCallback.h"
#include "composite_cpp.h"
#include "composite_cpp_terminate.h"
#include "gencoswin.h"
#include "rt_nonfinite.h"
#include <AudioFile.h> 


static void main_composite_cpp()
{

  // load clean speech
  AudioFile<float> clean_speech;

  clean_speech.load("../../data/clear8.wav");

  // load enhanced speech
  AudioFile<float> enhanced_speech;

  enhanced_speech.load("../../data/enhanced8.wav");


  t2r::Composite comp;

  float *sp_in = new float[clean_speech.getNumSamplesPerChannel()];
  float *enh_in = new float[enhanced_speech.getNumSamplesPerChannel()];

  for(int i = 0; i < clean_speech.getNumSamplesPerChannel(); i++){
    sp_in[i] = clean_speech.samples[0][i];
  }

  for(int i = 0; i < enhanced_speech.getNumSamplesPerChannel(); i++){
    enh_in[i] = enhanced_speech.samples[0][i];
  }

  t2r::OutputStructure out = comp(sp_in, clean_speech.getNumSamplesPerChannel(), (double)clean_speech.getSampleRate(), enh_in, enhanced_speech.getNumSamplesPerChannel(), (double)enhanced_speech.getSampleRate());

  std::cout << out.csig << " | " << out.cbak << " | " << out.covl << std::endl;

  delete sp_in;
  delete enh_in;

}

int main(int, const char * const [])
{
  // The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. 
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_composite_cpp();

  return 0;
}

// End of code generation (main.cpp)
