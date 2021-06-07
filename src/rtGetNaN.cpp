//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  rtGetNaN.cpp
//
//  Code generation for function 'composite_cpp'
//


//
//  Abstract:
//        MATLAB for code generation function to initialize non-finite, NaN

#include <rtGetNaN.h>

// Function: rtGetNaN ======================================================================
//  Abstract:
//  Initialize rtNaN needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE.

real_T rtGetNaN(void)
{
  return rtNaN;
}

// Function: rtGetNaNF =====================================================================
//  Abstract:
//  Initialize rtNaNF needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE.

real32_T rtGetNaNF(void)
{
  return rtNaNF;
}

// End of code generation (rtGetNaN.cpp)
