/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xgeqp3.c
 *
 * Code generation for function 'xgeqp3'
 *
 */

/* Include files */
#include "xgeqp3.h"
#include "only_bf_rtwutil.h"
#include "only_bf_types.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
int xgeqp3(emxArray_real_T *A, double tau_data[], int jpvt[12])
{
  double vn1[12];
  double vn2[12];
  double work[12];
  double d;
  double *A_data;
  int i;
  int k;
  int knt;
  int m;
  int ma;
  int pvt;
  int tau_size;
  A_data = A->data;
  m = A->size[0] - 1;
  tau_size = 12;
  ma = A->size[0];
  for (k = 0; k < 12; k++) {
    tau_data[k] = 0.0;
    jpvt[k] = k + 1;
    work[k] = 0.0;
    d = xnrm2(m + 1, A, k * ma + 1);
    vn1[k] = d;
    vn2[k] = d;
  }
  for (i = 0; i < 12; i++) {
    double atmp;
    double s;
    double smax;
    int b_i;
    int ii;
    int ip1;
    int ix;
    int lastc;
    int mmi;
    ip1 = i + 2;
    lastc = i * ma;
    ii = lastc + i;
    mmi = m - i;
    ix = 12 - i;
    knt = 0;
    if (12 - i > 1) {
      smax = fabs(vn1[i]);
      for (k = 2; k <= ix; k++) {
        s = fabs(vn1[(i + k) - 1]);
        if (s > smax) {
          knt = k - 1;
          smax = s;
        }
      }
    }
    pvt = i + knt;
    if (pvt != i) {
      ix = pvt * ma;
      for (k = 0; k <= m; k++) {
        knt = ix + k;
        smax = A_data[knt];
        b_i = lastc + k;
        A_data[knt] = A_data[b_i];
        A_data[b_i] = smax;
      }
      ix = jpvt[pvt];
      jpvt[pvt] = jpvt[i];
      jpvt[i] = ix;
      vn1[pvt] = vn1[i];
      vn2[pvt] = vn2[i];
    }
    atmp = A_data[ii];
    ix = ii + 2;
    tau_data[i] = 0.0;
    smax = xnrm2(mmi, A, ii + 2);
    if (smax != 0.0) {
      s = rt_hypotd_snf(A_data[ii], smax);
      if (A_data[ii] >= 0.0) {
        s = -s;
      }
      if (fabs(s) < 1.0020841800044864E-292) {
        knt = 0;
        b_i = ii + mmi;
        do {
          knt++;
          for (k = ix; k <= b_i + 1; k++) {
            A_data[k - 1] *= 9.9792015476736E+291;
          }
          s *= 9.9792015476736E+291;
          atmp *= 9.9792015476736E+291;
        } while ((fabs(s) < 1.0020841800044864E-292) && (knt < 20));
        s = rt_hypotd_snf(atmp, xnrm2(mmi, A, ii + 2));
        if (atmp >= 0.0) {
          s = -s;
        }
        tau_data[i] = (s - atmp) / s;
        smax = 1.0 / (atmp - s);
        for (k = ix; k <= b_i + 1; k++) {
          A_data[k - 1] *= smax;
        }
        for (k = 0; k < knt; k++) {
          s *= 1.0020841800044864E-292;
        }
        atmp = s;
      } else {
        tau_data[i] = (s - A_data[ii]) / s;
        smax = 1.0 / (A_data[ii] - s);
        b_i = ii + mmi;
        for (k = ix; k <= b_i + 1; k++) {
          A_data[k - 1] *= smax;
        }
        atmp = s;
      }
    }
    A_data[ii] = atmp;
    if (i + 1 < 12) {
      int jA;
      int lastv;
      atmp = A_data[ii];
      A_data[ii] = 1.0;
      jA = (ii + ma) + 1;
      if (tau_data[i] != 0.0) {
        bool exitg2;
        lastv = mmi;
        ix = ii + mmi;
        while ((lastv + 1 > 0) && (A_data[ix] == 0.0)) {
          lastv--;
          ix--;
        }
        lastc = 11 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          int exitg1;
          ix = jA + (lastc - 1) * ma;
          k = ix;
          do {
            exitg1 = 0;
            if (k <= ix + lastv) {
              if (A_data[k - 1] != 0.0) {
                exitg1 = 1;
              } else {
                k++;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
        lastc--;
      } else {
        lastv = -1;
        lastc = -1;
      }
      if (lastv + 1 > 0) {
        if (lastc + 1 != 0) {
          memset(&work[0], 0, (unsigned int)(lastc + 1) * sizeof(double));
          knt = 0;
          b_i = jA + ma * lastc;
          for (pvt = jA; ma < 0 ? pvt >= b_i : pvt <= b_i; pvt += ma) {
            smax = 0.0;
            ix = pvt + lastv;
            for (k = pvt; k <= ix; k++) {
              smax += A_data[k - 1] * A_data[(ii + k) - pvt];
            }
            work[knt] += smax;
            knt++;
          }
        }
        if (!(-tau_data[i] == 0.0)) {
          for (pvt = 0; pvt <= lastc; pvt++) {
            d = work[pvt];
            if (d != 0.0) {
              smax = d * -tau_data[i];
              b_i = lastv + jA;
              for (knt = jA; knt <= b_i; knt++) {
                A_data[knt - 1] += A_data[(ii + knt) - jA] * smax;
              }
            }
            jA += ma;
          }
        }
      }
      A_data[ii] = atmp;
    }
    for (pvt = ip1; pvt < 13; pvt++) {
      ix = i + (pvt - 1) * ma;
      d = vn1[pvt - 1];
      if (d != 0.0) {
        smax = fabs(A_data[ix]) / d;
        smax = 1.0 - smax * smax;
        if (smax < 0.0) {
          smax = 0.0;
        }
        s = d / vn2[pvt - 1];
        s = smax * (s * s);
        if (s <= 1.4901161193847656E-8) {
          d = xnrm2(mmi, A, ix + 2);
          vn1[pvt - 1] = d;
          vn2[pvt - 1] = d;
        } else {
          vn1[pvt - 1] = d * sqrt(smax);
        }
      }
    }
  }
  return tau_size;
}

/* End of code generation (xgeqp3.c) */
