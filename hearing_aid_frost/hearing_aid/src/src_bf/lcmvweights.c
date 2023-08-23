/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lcmvweights.c
 *
 * Code generation for function 'lcmvweights'
 *
 */

/* Include files */
#include "lcmvweights.h"
#include "only_bf_emxutil.h"
#include "only_bf_types.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "xgeqp3.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void qrlinsolve(const emxArray_real_T *A, double Xout[24], double Fout[24])
{
  static const signed char B[24] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
  emxArray_real_T *b_A;
  double R_data[144];
  double c_A_data[144];
  double X[24];
  double X_data[24];
  double tau_data[12];
  const double *A_data;
  double s;
  double smax;
  double *b_A_data;
  int jpvt[12];
  int a;
  int b_i;
  int b_tmp;
  int i;
  int i1;
  int j;
  int jA;
  int jp1j;
  int k;
  int kAcol;
  int mmj_tmp;
  A_data = A->data;
  emxInit_real_T(&b_A);
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[1];
  b_A->size[1] = 12;
  emxEnsureCapacity_real_T(b_A, i);
  b_A_data = b_A->data;
  jA = A->size[1];
  for (i = 0; i < 12; i++) {
    for (i1 = 0; i1 < jA; i1++) {
      b_A_data[i1 + b_A->size[0] * i] = A_data[i + 12 * i1];
    }
  }
  xgeqp3(b_A, tau_data, jpvt);
  b_A_data = b_A->data;
  for (j = 0; j < 12; j++) {
    for (b_i = 0; b_i <= j; b_i++) {
      R_data[b_i + 12 * j] = b_A_data[b_i + b_A->size[0] * j];
    }
    i = j + 2;
    if (i <= 12) {
      memset(&R_data[(j * 12 + i) + -1], 0,
             (unsigned int)(-i + 13) * sizeof(double));
    }
    tau_data[j] = jpvt[j];
  }
  emxFree_real_T(&b_A);
  for (i = 0; i < 12; i++) {
    for (i1 = 0; i1 < 12; i1++) {
      c_A_data[i1 + 12 * i] = R_data[i + 12 * i1];
    }
    jA = (int)tau_data[i];
    X[i] = B[jA - 1];
    X[i + 12] = B[jA + 11];
    jpvt[i] = i + 1;
  }
  for (j = 0; j < 11; j++) {
    mmj_tmp = 10 - j;
    b_tmp = j * 13;
    jp1j = b_tmp + 2;
    jA = 12 - j;
    a = 0;
    smax = fabs(c_A_data[b_tmp]);
    for (k = 2; k <= jA; k++) {
      s = fabs(c_A_data[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (c_A_data[b_tmp + a] != 0.0) {
      if (a != 0) {
        jA = j + a;
        jpvt[j] = jA + 1;
        for (k = 0; k < 12; k++) {
          a = j + k * 12;
          smax = c_A_data[a];
          i = jA + k * 12;
          c_A_data[a] = c_A_data[i];
          c_A_data[i] = smax;
        }
      }
      i = (b_tmp - j) + 12;
      for (b_i = jp1j; b_i <= i; b_i++) {
        c_A_data[b_i - 1] /= c_A_data[b_tmp];
      }
    }
    jA = b_tmp;
    for (a = 0; a <= mmj_tmp; a++) {
      smax = c_A_data[(b_tmp + a * 12) + 12];
      if (smax != 0.0) {
        i = jA + 14;
        i1 = (jA - j) + 24;
        for (kAcol = i; kAcol <= i1; kAcol++) {
          c_A_data[kAcol - 1] += c_A_data[((b_tmp + kAcol) - jA) - 13] * -smax;
        }
      }
      jA += 12;
    }
    i = jpvt[j];
    if (i != j + 1) {
      jA = (int)X[j];
      X[j] = X[i - 1];
      X[i - 1] = jA;
      jA = (int)X[j + 12];
      X[j + 12] = X[i + 11];
      X[i + 11] = jA;
    }
  }
  for (j = 0; j < 2; j++) {
    jp1j = 12 * j;
    for (k = 0; k < 12; k++) {
      kAcol = 12 * k;
      i = k + jp1j;
      if (X[i] != 0.0) {
        i1 = k + 2;
        for (b_i = i1; b_i < 13; b_i++) {
          jA = (b_i + jp1j) - 1;
          X[jA] -= X[i] * c_A_data[(b_i + kAcol) - 1];
        }
      }
    }
  }
  for (j = 0; j < 2; j++) {
    jp1j = 12 * j;
    for (k = 11; k >= 0; k--) {
      kAcol = 12 * k;
      i = k + jp1j;
      smax = X[i];
      if (smax != 0.0) {
        X[i] = smax / c_A_data[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          jA = b_i + jp1j;
          X[jA] -= X[i] * c_A_data[b_i + kAcol];
        }
      }
    }
  }
  memcpy(&X_data[0], &X[0], 24U * sizeof(double));
  for (i = 0; i < 12; i++) {
    jpvt[i] = i + 1;
  }
  for (j = 0; j < 11; j++) {
    mmj_tmp = 10 - j;
    b_tmp = j * 13;
    jp1j = b_tmp + 2;
    jA = 12 - j;
    a = 0;
    smax = fabs(R_data[b_tmp]);
    for (k = 2; k <= jA; k++) {
      s = fabs(R_data[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (R_data[b_tmp + a] != 0.0) {
      if (a != 0) {
        jA = j + a;
        jpvt[j] = jA + 1;
        for (k = 0; k < 12; k++) {
          a = j + k * 12;
          smax = R_data[a];
          i = jA + k * 12;
          R_data[a] = R_data[i];
          R_data[i] = smax;
        }
      }
      i = (b_tmp - j) + 12;
      for (b_i = jp1j; b_i <= i; b_i++) {
        R_data[b_i - 1] /= R_data[b_tmp];
      }
    }
    jA = b_tmp;
    for (a = 0; a <= mmj_tmp; a++) {
      smax = R_data[(b_tmp + a * 12) + 12];
      if (smax != 0.0) {
        i = jA + 14;
        i1 = (jA - j) + 24;
        for (kAcol = i; kAcol <= i1; kAcol++) {
          R_data[kAcol - 1] += R_data[((b_tmp + kAcol) - jA) - 13] * -smax;
        }
      }
      jA += 12;
    }
    i = jpvt[j];
    if (i != j + 1) {
      smax = X_data[j];
      X_data[j] = X_data[i - 1];
      X_data[i - 1] = smax;
      smax = X_data[j + 12];
      X_data[j + 12] = X_data[i + 11];
      X_data[i + 11] = smax;
    }
  }
  for (j = 0; j < 2; j++) {
    jp1j = 12 * j;
    for (k = 0; k < 12; k++) {
      kAcol = 12 * k;
      i = k + jp1j;
      if (X_data[i] != 0.0) {
        i1 = k + 2;
        for (b_i = i1; b_i < 13; b_i++) {
          a = (b_i + jp1j) - 1;
          X_data[a] -= X_data[i] * R_data[(b_i + kAcol) - 1];
        }
      }
    }
  }
  for (j = 0; j < 2; j++) {
    jp1j = 12 * j;
    for (k = 11; k >= 0; k--) {
      kAcol = 12 * k;
      i = k + jp1j;
      smax = X_data[i];
      if (smax != 0.0) {
        X_data[i] = smax / R_data[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          i1 = b_i + jp1j;
          X_data[i1] -= X_data[i] * R_data[b_i + kAcol];
        }
      }
    }
  }
  sort(tau_data, jpvt);
  for (i = 0; i < 12; i++) {
    i1 = jpvt[i];
    Xout[i] = X_data[i1 - 1];
    Fout[i] = X[i1 - 1];
    Xout[i + 12] = X_data[i1 + 11];
    Fout[i + 12] = X[i1 + 11];
  }
}

/* End of code generation (lcmvweights.c) */
