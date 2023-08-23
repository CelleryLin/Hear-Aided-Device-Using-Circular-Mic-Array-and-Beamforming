/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * only_bf_rtwutil.c
 *
 * Code generation for function 'only_bf_rtwutil'
 *
 */

/* Include files */
#include "only_bf_rtwutil.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

/* End of code generation (only_bf_rtwutil.c) */
