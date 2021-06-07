//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  relop.cpp
//
//  Code generation for function 'relop'
//


// Include files
#include <relop.h>
#include <composite_cpp.h>
#include <rt_nonfinite.h>
#include <cmath>
#include <math.h>

// Function Definitions
boolean_T iseq(double x, double y)
{
  boolean_T p;
  double absx;
  int exponent;
  absx = std::abs(y / 2.0);
  if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = std::ldexp(1.0, exponent - 53);
    }
  } else {
    absx = rtNaN;
  }

  if ((std::abs(y - x) < absx) || (rtIsInf(x) && rtIsInf(y) && ((x > 0.0) == (y >
         0.0)))) {
    p = true;
  } else {
    p = false;
  }

  return p;
}

// End of code generation (relop.cpp)
