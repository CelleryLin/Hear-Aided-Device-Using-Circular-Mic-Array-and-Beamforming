/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * BlockSMIWeightsEstimator.c
 *
 * Code generation for function 'BlockSMIWeightsEstimator'
 *
 */

/* Include files */
#include "BlockSMIWeightsEstimator.h"
#include "lcmvweights.h"
#include "only_bf_emxutil.h"
#include "only_bf_internal_types.h"
#include "only_bf_rtwutil.h"
#include "only_bf_types.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "xnrm2.h"
#include <math.h>

/* Function Declarations */
static int div_nde_s32_floor(int numerator);

/* Function Definitions */
static int div_nde_s32_floor(int numerator)
{
  int i;
  if ((numerator < 0) && (numerator % 12 != 0)) {
    i = -1;
  } else {
    i = 0;
  }
  return numerator / 12 + i;
}

void c_BlockSMIWeightsEstimator_step(const c_phased_internal_BlockSMIWeigh *obj,
                                     const double x[14400], double w[12])
{
  static const signed char b[144] = {
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  emxArray_real_T *b_x;
  emxArray_real_T *c_x;
  double xArg[14400];
  double F[24];
  double temp[24];
  double A[4];
  double R[4];
  double tau[2];
  double vn1[2];
  double vn2[2];
  double work[2];
  double absxk;
  double delta;
  double scale;
  double t;
  double *b_x_data;
  double *x_data;
  int jpvt[2];
  int b_i;
  int c_i;
  int i;
  int iy;
  int k;
  int kend;
  int pvt;
  for (i = 0; i < 14400; i++) {
    xArg[i] = x[i] / 34.641016151377549;
  }
  delta = sqrt(obj->DiagonalLoadingFactor);
  emxInit_real_T(&b_x);
  if (delta != 0.0) {
    i = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1212;
    b_x->size[1] = 12;
    emxEnsureCapacity_real_T(b_x, i);
    x_data = b_x->data;
    for (i = 0; i < 12; i++) {
      for (pvt = 0; pvt < 1200; pvt++) {
        x_data[pvt + b_x->size[0] * i] = xArg[pvt + 1200 * i];
      }
      for (pvt = 0; pvt < 12; pvt++) {
        x_data[(pvt + b_x->size[0] * i) + 1200] =
            delta * (double)b[pvt + 12 * i];
      }
    }
  } else {
    i = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1200;
    b_x->size[1] = 12;
    emxEnsureCapacity_real_T(b_x, i);
    x_data = b_x->data;
    for (i = 0; i < 14400; i++) {
      x_data[i] = xArg[i];
    }
  }
  emxInit_real_T(&c_x);
  i = c_x->size[0] * c_x->size[1];
  c_x->size[0] = 12;
  c_x->size[1] = b_x->size[0];
  emxEnsureCapacity_real_T(c_x, i);
  b_x_data = c_x->data;
  kend = b_x->size[0];
  for (i = 0; i < kend; i++) {
    for (pvt = 0; pvt < 12; pvt++) {
      b_x_data[pvt + 12 * i] = x_data[i + b_x->size[0] * pvt];
    }
  }
  emxFree_real_T(&b_x);
  qrlinsolve(c_x, temp, F);
  emxFree_real_T(&c_x);
  for (k = 0; k < 2; k++) {
    jpvt[k] = k + 1;
    tau[k] = 0.0;
    work[k] = 0.0;
    iy = k * 12;
    delta = 0.0;
    scale = 3.3121686421112381E-170;
    kend = iy + 12;
    for (c_i = iy + 1; c_i <= kend; c_i++) {
      absxk = fabs(F[c_i - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        delta = delta * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        delta += t * t;
      }
    }
    t = scale * sqrt(delta);
    vn1[k] = t;
    vn2[k] = t;
  }
  for (b_i = 0; b_i < 2; b_i++) {
    int ii;
    int ip1;
    int lastv;
    ip1 = b_i + 2;
    ii = b_i * 12 + b_i;
    iy = 0;
    if ((2 - b_i > 1) && (fabs(vn1[b_i + 1]) > fabs(vn1[b_i]))) {
      iy = 1;
    }
    pvt = b_i + iy;
    if (pvt != b_i) {
      kend = pvt * 12;
      iy = b_i * 12;
      for (k = 0; k < 12; k++) {
        c_i = kend + k;
        delta = F[c_i];
        lastv = iy + k;
        F[c_i] = F[lastv];
        F[lastv] = delta;
      }
      iy = jpvt[pvt];
      jpvt[pvt] = jpvt[b_i];
      jpvt[b_i] = iy;
      vn1[pvt] = vn1[b_i];
      vn2[pvt] = vn2[b_i];
    }
    absxk = F[ii];
    iy = ii + 2;
    tau[b_i] = 0.0;
    delta = b_xnrm2(11 - b_i, F, ii + 2);
    if (delta != 0.0) {
      t = F[ii];
      scale = rt_hypotd_snf(t, delta);
      if (t >= 0.0) {
        scale = -scale;
      }
      if (fabs(scale) < 1.0020841800044864E-292) {
        kend = 0;
        i = (ii - b_i) + 12;
        do {
          kend++;
          for (k = iy; k <= i; k++) {
            F[k - 1] *= 9.9792015476736E+291;
          }
          scale *= 9.9792015476736E+291;
          absxk *= 9.9792015476736E+291;
        } while ((fabs(scale) < 1.0020841800044864E-292) && (kend < 20));
        scale = rt_hypotd_snf(absxk, b_xnrm2(11 - b_i, F, ii + 2));
        if (absxk >= 0.0) {
          scale = -scale;
        }
        tau[b_i] = (scale - absxk) / scale;
        delta = 1.0 / (absxk - scale);
        for (k = iy; k <= i; k++) {
          F[k - 1] *= delta;
        }
        for (k = 0; k < kend; k++) {
          scale *= 1.0020841800044864E-292;
        }
        absxk = scale;
      } else {
        tau[b_i] = (scale - t) / scale;
        delta = 1.0 / (t - scale);
        i = (ii - b_i) + 12;
        for (k = iy; k <= i; k++) {
          F[k - 1] *= delta;
        }
        absxk = scale;
      }
    }
    F[ii] = absxk;
    if (b_i + 1 < 2) {
      F[ii] = 1.0;
      if (tau[0] != 0.0) {
        lastv = 12;
        c_i = ii + 11;
        while ((lastv > 0) && (F[c_i] == 0.0)) {
          lastv--;
          c_i--;
        }
        iy = 1;
        kend = ii + 12;
        int exitg1;
        do {
          exitg1 = 0;
          if (kend + 1 <= (ii + lastv) + 12) {
            if (F[kend] != 0.0) {
              exitg1 = 1;
            } else {
              kend++;
            }
          } else {
            iy = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      } else {
        lastv = 0;
        iy = 0;
      }
      if (lastv > 0) {
        c_i = ii + 13;
        if (iy != 0) {
          work[0] = 0.0;
          for (pvt = c_i; pvt <= c_i; pvt += 12) {
            delta = 0.0;
            i = (pvt + lastv) - 1;
            for (kend = pvt; kend <= i; kend++) {
              delta += F[kend - 1] * F[(ii + kend) - pvt];
            }
            kend = div_nde_s32_floor((pvt - ii) - 13);
            work[kend] += delta;
          }
        }
        if (!(-tau[0] == 0.0)) {
          kend = ii;
          for (k = 0; k < iy; k++) {
            if (work[0] != 0.0) {
              delta = work[0] * -tau[0];
              i = kend + 13;
              pvt = lastv + kend;
              for (c_i = i; c_i <= pvt + 12; c_i++) {
                F[c_i - 1] += F[((ii + c_i) - kend) - 13] * delta;
              }
            }
            kend += 12;
          }
        }
      }
      F[ii] = absxk;
    }
    for (k = ip1; k < 3; k++) {
      if (vn1[1] != 0.0) {
        delta = fabs(F[b_i + 12]) / vn1[1];
        delta = 1.0 - delta * delta;
        if (delta < 0.0) {
          delta = 0.0;
        }
        scale = vn1[1] / vn2[1];
        scale = delta * (scale * scale);
        if (scale <= 1.4901161193847656E-8) {
          t = b_xnrm2(11 - b_i, F, b_i + 14);
          vn1[1] = t;
          vn2[1] = t;
        } else {
          vn1[1] *= sqrt(delta);
        }
      }
    }
  }
  for (k = 0; k < 2; k++) {
    for (b_i = 0; b_i <= k; b_i++) {
      R[b_i + (k << 1)] = F[b_i + 12 * k];
    }
    if (k + 2 <= 2) {
      R[(k << 1) + 1] = 0.0;
    }
    kend = jpvt[k];
    work[k] = kend;
    tau[k] = kend;
  }
  b_sort(tau, jpvt);
  A[0] = R[0];
  A[1] = R[2];
  A[2] = R[1];
  A[3] = R[3];
  t = fabs(R[0]);
  if (fabs(R[2]) > t) {
    kend = 1;
    iy = 0;
  } else {
    kend = 0;
    iy = 1;
  }
  scale = A[iy] / A[kend];
  c_i = 2 - (int)work[kend];
  delta = A[kend + 2];
  vn1[1] =
      ((2.0 - work[iy]) - (double)c_i * scale) / (A[iy + 2] - scale * delta);
  vn1[0] = ((double)c_i - vn1[1] * delta) / A[kend];
  if (fabs(R[1]) > t) {
    kend = 1;
    iy = 0;
  } else {
    kend = 0;
    iy = 1;
  }
  scale = R[iy] / R[kend];
  delta = R[kend + 2];
  tau[1] = (vn1[iy] - vn1[kend] * scale) / (R[iy + 2] - scale * delta);
  tau[0] = (vn1[kend] - tau[1] * delta) / R[kend];
  t = tau[jpvt[0] - 1];
  delta = tau[jpvt[1] - 1];
  for (i = 0; i < 12; i++) {
    w[i] = temp[i] * t + temp[i + 12] * delta;
  }
}

/* End of code generation (BlockSMIWeightsEstimator.c) */
