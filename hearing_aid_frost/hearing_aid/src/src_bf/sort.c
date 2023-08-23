/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void b_sort(double x[2], int idx[2])
{
  if ((x[0] <= x[1]) || rtIsNaN(x[1])) {
    idx[0] = 1;
    idx[1] = 2;
  } else {
    double tmp;
    idx[0] = 2;
    idx[1] = 1;
    tmp = x[0];
    x[0] = x[1];
    x[1] = tmp;
  }
}

void sort(double x[12], int idx[12])
{
  double xwork[12];
  double x4[4];
  int i1;
  int i3;
  int ib;
  int k;
  int nNaNs;
  int quartetOffset;
  signed char idx4[4];
  x4[0] = 0.0;
  idx4[0] = 0;
  x4[1] = 0.0;
  idx4[1] = 0;
  x4[2] = 0.0;
  idx4[2] = 0;
  x4[3] = 0.0;
  idx4[3] = 0;
  memset(&xwork[0], 0, 12U * sizeof(double));
  for (ib = 0; ib < 12; ib++) {
    idx[ib] = 0;
  }
  nNaNs = 0;
  ib = 0;
  for (k = 0; k < 12; k++) {
    if (rtIsNaN(x[k])) {
      idx[11 - nNaNs] = k + 1;
      xwork[11 - nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (signed char)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        double d;
        double d1;
        int i4;
        signed char b_i1;
        signed char b_i3;
        signed char i;
        signed char i2;
        quartetOffset = k - nNaNs;
        if (x4[0] <= x4[1]) {
          i1 = 1;
          ib = 2;
        } else {
          i1 = 2;
          ib = 1;
        }
        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }
        d = x4[i3 - 1];
        d1 = x4[i1 - 1];
        if (d1 <= d) {
          d1 = x4[ib - 1];
          if (d1 <= d) {
            i = (signed char)i1;
            b_i1 = (signed char)ib;
            i2 = (signed char)i3;
            b_i3 = (signed char)i4;
          } else if (d1 <= x4[i4 - 1]) {
            i = (signed char)i1;
            b_i1 = (signed char)i3;
            i2 = (signed char)ib;
            b_i3 = (signed char)i4;
          } else {
            i = (signed char)i1;
            b_i1 = (signed char)i3;
            i2 = (signed char)i4;
            b_i3 = (signed char)ib;
          }
        } else {
          d = x4[i4 - 1];
          if (d1 <= d) {
            if (x4[ib - 1] <= d) {
              i = (signed char)i3;
              b_i1 = (signed char)i1;
              i2 = (signed char)ib;
              b_i3 = (signed char)i4;
            } else {
              i = (signed char)i3;
              b_i1 = (signed char)i1;
              i2 = (signed char)i4;
              b_i3 = (signed char)ib;
            }
          } else {
            i = (signed char)i3;
            b_i1 = (signed char)i4;
            i2 = (signed char)i1;
            b_i3 = (signed char)ib;
          }
        }
        idx[quartetOffset - 3] = idx4[i - 1];
        idx[quartetOffset - 2] = idx4[b_i1 - 1];
        idx[quartetOffset - 1] = idx4[i2 - 1];
        idx[quartetOffset] = idx4[b_i3 - 1];
        x[quartetOffset - 3] = x4[i - 1];
        x[quartetOffset - 2] = x4[b_i1 - 1];
        x[quartetOffset - 1] = x4[i2 - 1];
        x[quartetOffset] = x4[b_i3 - 1];
        ib = 0;
      }
    }
  }
  if (ib > 0) {
    signed char perm[4];
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }
    quartetOffset = (unsigned char)ib;
    for (k = 0; k < quartetOffset; k++) {
      i1 = perm[k] - 1;
      i3 = ((k - nNaNs) - ib) + 12;
      idx[i3] = idx4[i1];
      x[i3] = x4[i1];
    }
  }
  ib = (nNaNs >> 1) + 12;
  for (k = 0; k <= ib - 13; k++) {
    i1 = (k - nNaNs) + 12;
    quartetOffset = idx[i1];
    idx[i1] = idx[11 - k];
    idx[11 - k] = quartetOffset;
    x[i1] = xwork[11 - k];
    x[11 - k] = xwork[i1];
  }
  if ((nNaNs & 1) != 0) {
    quartetOffset = ib - nNaNs;
    x[quartetOffset] = xwork[quartetOffset];
  }
  if (12 - nNaNs > 1) {
    int iwork[12];
    for (ib = 0; ib < 12; ib++) {
      iwork[ib] = 0;
    }
    i3 = (12 - nNaNs) >> 2;
    quartetOffset = 4;
    while (i3 > 1) {
      if ((i3 & 1) != 0) {
        i3--;
        ib = quartetOffset * i3;
        i1 = 12 - (nNaNs + ib);
        if (i1 > quartetOffset) {
          merge(idx, x, ib, quartetOffset, i1 - quartetOffset, iwork, xwork);
        }
      }
      i3 >>= 1;
      for (k = 0; k < i3; k++) {
        merge(idx, x, 0, quartetOffset, quartetOffset, iwork, xwork);
      }
      quartetOffset <<= 1;
    }
    if (12 - nNaNs > quartetOffset) {
      merge(idx, x, 0, quartetOffset, 12 - (nNaNs + quartetOffset), iwork,
            xwork);
    }
  }
}

/* End of code generation (sort.c) */
