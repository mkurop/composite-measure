//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  main.cpp
//
//  Code generation for function 'main'
//


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
#include "main.h"
#include "FFTImplementationCallback.h"
#include "composite_cpp.h"
#include "composite_cpp_terminate.h"
#include "gencoswin.h"
#include "rt_nonfinite.h"

// Function Declarations
static coder::array<double, 1U> argInit_Unboundedx1_real_T();
static double argInit_real_T();
static void main_composite_cpp();

// Function Definitions
static coder::array<double, 1U> argInit_Unboundedx1_real_T()
{
  coder::array<double, 1U> result;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result.set_size(2);

  // Loop over the array to initialize each element.
  for (int idx0 = 0; idx0 < result.size(0); idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx0] = argInit_real_T();
  }

  return result;
}

static double argInit_real_T()
{
  return 0.0;
}

static void main_composite_cpp()
{
  coder::array<double, 1U> cleanFile_tmp;
  double Srate1_tmp;
  double Csig_data[2];
  int Csig_size[2];
  double Cbak_data[2];
  int Cbak_size[2];
  double Covl_data[2];
  int Covl_size[2];

  // Initialize function 'composite_cpp' input arguments.
  // Initialize function input argument 'cleanFile'.
  cleanFile_tmp = argInit_Unboundedx1_real_T();

  // Initialize function input argument 'enhancedFile'.
  Srate1_tmp = argInit_real_T();

  // Call the entry-point 'composite_cpp'.
  composite_cpp(cleanFile_tmp, cleanFile_tmp, Srate1_tmp, Srate1_tmp, Csig_data,
                Csig_size, Cbak_data, Cbak_size, Covl_data, Covl_size);
}

int main(int, const char * const [])
{
  // The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here. 
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_composite_cpp();

  // Terminate the application.
  // You do not need to do this more than one time.
  composite_cpp_terminate();
  return 0;
}

// End of code generation (main.cpp)
