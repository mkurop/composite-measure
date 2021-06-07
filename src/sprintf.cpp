//
//  Trial License - for use to evaluate programs for possible purchase as
//  an end-user only.
//
//  sprintf.cpp
//
//  Code generation for function 'sprintf'
//


// Include files
#include <sprintf.h>
#include <composite_cpp.h>
#include <rt_nonfinite.h>
#include <stdio.h>

// Type Definitions
struct emxArray_char_T_1x11
{
  char data[11];
  int size[2];
};

struct cell_wrap_6
{
  emxArray_char_T_1x11 f1;
};

// Function Definitions
void b_sprintf(const char varargin_1_data[], const int varargin_1_size[2], coder::
               array<char, 2U> &str)
{
  cell_wrap_6 validatedHoleFilling[1];
  int nbytes;
  int i;
  char b_varargin_1_data[11];
  char c_varargin_1_data[11];
  validatedHoleFilling[0].f1.size[1] = varargin_1_size[1] + 1;
  nbytes = varargin_1_size[1];
  for (i = 0; i < nbytes; i++) {
    validatedHoleFilling[0].f1.data[i] = varargin_1_data[i];
  }

  validatedHoleFilling[0].f1.data[varargin_1_size[1]] = '\x00';
  nbytes = validatedHoleFilling[0].f1.size[1];
  for (i = 0; i < nbytes; i++) {
    b_varargin_1_data[i] = validatedHoleFilling[0].f1.data[i];
  }

  nbytes = validatedHoleFilling[0].f1.size[1];
  for (i = 0; i < nbytes; i++) {
    c_varargin_1_data[i] = validatedHoleFilling[0].f1.data[i];
  }

  nbytes = snprintf(NULL, 0, "Mode: \"%s\" is unsupported.", &c_varargin_1_data
                    [0]);
  str.set_size(1, (nbytes + 1));
  snprintf(&str[0], (size_t)(nbytes + 1), "Mode: \"%s\" is unsupported.",
           &b_varargin_1_data[0]);
  if (1 > nbytes) {
    nbytes = 0;
  }

  str.set_size(str.size(0), nbytes);
}

// End of code generation (sprintf.cpp)
