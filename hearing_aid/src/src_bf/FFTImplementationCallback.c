/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * FFTImplementationCallback.c
 *
 * Code generation for function 'FFTImplementationCallback'
 *
 */

/* Include files */
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void c_FFTImplementationCallback_doH(const double x[1200], creal_T y_data[],
                                     const double costab_data[],
                                     const double sintab_data[])
{
  creal_T reconVar1_data[1024];
  creal_T reconVar2_data[1024];
  double hcostab_data[512];
  double hsintab_data[512];
  double im;
  double re;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_im_tmp;
  double temp_re;
  double temp_re_tmp;
  int bitrevIndex_data[1024];
  int wrapIndex_data[1024];
  int i;
  int iDelta;
  int iDelta2;
  int iheight;
  int ihi;
  int iy;
  int ju;
  int k;
  for (i = 0; i < 512; i++) {
    iy = ((i + 1) << 1) - 2;
    hcostab_data[i] = costab_data[iy];
    hsintab_data[i] = sintab_data[iy];
  }
  ju = 0;
  iy = 1;
  for (i = 0; i < 1024; i++) {
    temp2_re = sintab_data[i];
    reconVar1_data[i].re = temp2_re + 1.0;
    temp2_im = costab_data[i];
    reconVar1_data[i].im = -temp2_im;
    reconVar2_data[i].re = 1.0 - temp2_re;
    reconVar2_data[i].im = temp2_im;
    if (i + 1 != 1) {
      wrapIndex_data[i] = 1025 - i;
    } else {
      wrapIndex_data[0] = 1;
    }
    bitrevIndex_data[i] = 0;
  }
  for (ihi = 0; ihi < 1023; ihi++) {
    bool tst;
    bitrevIndex_data[ihi] = iy;
    iy = 1024;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju + 1;
  }
  bitrevIndex_data[1023] = iy;
  for (i = 0; i < 600; i++) {
    iy = i << 1;
    ju = bitrevIndex_data[i];
    y_data[ju - 1].re = x[iy];
    y_data[ju - 1].im = x[iy + 1];
  }
  for (i = 0; i <= 1022; i += 2) {
    temp2_re = y_data[i + 1].re;
    temp2_im = y_data[i + 1].im;
    temp_re = temp2_re;
    temp_im = temp2_im;
    re = y_data[i].re;
    im = y_data[i].im;
    temp2_re = re - temp2_re;
    temp2_im = im - temp2_im;
    y_data[i + 1].re = temp2_re;
    y_data[i + 1].im = temp2_im;
    re += temp_re;
    im += temp_im;
    y_data[i].re = re;
    y_data[i].im = im;
  }
  iDelta = 2;
  iDelta2 = 4;
  k = 256;
  iheight = 1021;
  while (k > 0) {
    for (i = 0; i < iheight; i += iDelta2) {
      iy = i + iDelta;
      temp_re = y_data[iy].re;
      temp_im = y_data[iy].im;
      y_data[iy].re = y_data[i].re - temp_re;
      y_data[iy].im = y_data[i].im - temp_im;
      y_data[i].re += temp_re;
      y_data[i].im += temp_im;
    }
    iy = 1;
    for (ju = k; ju < 512; ju += k) {
      temp2_re = hcostab_data[ju];
      temp2_im = hsintab_data[ju];
      i = iy;
      ihi = iy + iheight;
      while (i < ihi) {
        int temp_re_tmp_tmp;
        temp_re_tmp_tmp = i + iDelta;
        temp_re_tmp = y_data[temp_re_tmp_tmp].im;
        temp_im = y_data[temp_re_tmp_tmp].re;
        temp_re = temp2_re * temp_im - temp2_im * temp_re_tmp;
        temp_im = temp2_re * temp_re_tmp + temp2_im * temp_im;
        y_data[temp_re_tmp_tmp].re = y_data[i].re - temp_re;
        y_data[temp_re_tmp_tmp].im = y_data[i].im - temp_im;
        y_data[i].re += temp_re;
        y_data[i].im += temp_im;
        i += iDelta2;
      }
      iy++;
    }
    k /= 2;
    iDelta = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iDelta;
  }
  temp2_re = y_data[0].re;
  temp_im_tmp = y_data[0].im;
  temp_im = temp2_re * reconVar1_data[0].re;
  re = temp2_re * reconVar1_data[0].im;
  temp_re = temp2_re * reconVar2_data[0].re;
  temp2_im = temp2_re * reconVar2_data[0].im;
  y_data[0].re = 0.5 * ((temp_im - temp_im_tmp * reconVar1_data[0].im) +
                        (temp_re - -temp_im_tmp * reconVar2_data[0].im));
  y_data[0].im = 0.5 * ((re + temp_im_tmp * reconVar1_data[0].re) +
                        (temp2_im + -temp_im_tmp * reconVar2_data[0].re));
  y_data[1024].re = 0.5 * ((temp_re - temp_im_tmp * reconVar2_data[0].im) +
                           (temp_im - -temp_im_tmp * reconVar1_data[0].im));
  y_data[1024].im = 0.5 * ((temp2_im + temp_im_tmp * reconVar2_data[0].re) +
                           (re + -temp_im_tmp * reconVar1_data[0].re));
  for (i = 0; i < 511; i++) {
    double temp2_im_tmp;
    temp_re_tmp = y_data[i + 1].re;
    temp_im_tmp = y_data[i + 1].im;
    iy = wrapIndex_data[i + 1];
    temp2_im = y_data[iy - 1].re;
    temp2_im_tmp = y_data[iy - 1].im;
    temp_im = reconVar1_data[i + 1].im;
    temp_re = reconVar1_data[i + 1].re;
    re = reconVar2_data[i + 1].im;
    im = reconVar2_data[i + 1].re;
    y_data[i + 1].re = 0.5 * ((temp_re_tmp * temp_re - temp_im_tmp * temp_im) +
                              (temp2_im * im - -temp2_im_tmp * re));
    y_data[i + 1].im = 0.5 * ((temp_re_tmp * temp_im + temp_im_tmp * temp_re) +
                              (temp2_im * re + -temp2_im_tmp * im));
    y_data[i + 1025].re =
        0.5 * ((temp_re_tmp * im - temp_im_tmp * re) +
               (temp2_im * temp_re - -temp2_im_tmp * temp_im));
    y_data[i + 1025].im =
        0.5 * ((temp_re_tmp * re + temp_im_tmp * im) +
               (temp2_im * temp_im + -temp2_im_tmp * temp_re));
    re = reconVar1_data[iy - 1].im;
    im = reconVar1_data[iy - 1].re;
    temp_im = reconVar2_data[iy - 1].im;
    temp2_re = reconVar2_data[iy - 1].re;
    y_data[iy - 1].re =
        0.5 * ((temp2_im * im - temp2_im_tmp * re) +
               (temp_re_tmp * temp2_re - -temp_im_tmp * temp_im));
    y_data[iy - 1].im =
        0.5 * ((temp2_im * re + temp2_im_tmp * im) +
               (temp_re_tmp * temp_im + -temp_im_tmp * temp2_re));
    y_data[iy + 1023].re =
        0.5 * ((temp2_im * temp2_re - temp2_im_tmp * temp_im) +
               (temp_re_tmp * im - -temp_im_tmp * re));
    y_data[iy + 1023].im =
        0.5 * ((temp2_im * temp_im + temp2_im_tmp * temp2_re) +
               (temp_re_tmp * re + -temp_im_tmp * im));
  }
  temp2_re = y_data[512].re;
  temp_im_tmp = y_data[512].im;
  temp_im = temp2_re * reconVar1_data[512].re;
  re = temp2_re * reconVar1_data[512].im;
  temp_re = temp2_re * reconVar2_data[512].re;
  temp2_im = temp2_re * reconVar2_data[512].im;
  y_data[512].re = 0.5 * ((temp_im - temp_im_tmp * reconVar1_data[512].im) +
                          (temp_re - -temp_im_tmp * reconVar2_data[512].im));
  y_data[512].im = 0.5 * ((re + temp_im_tmp * reconVar1_data[512].re) +
                          (temp2_im + -temp_im_tmp * reconVar2_data[512].re));
  y_data[1536].re = 0.5 * ((temp_re - temp_im_tmp * reconVar2_data[512].im) +
                           (temp_im - -temp_im_tmp * reconVar1_data[512].im));
  y_data[1536].im = 0.5 * ((temp2_im + temp_im_tmp * reconVar2_data[512].re) +
                           (re + -temp_im_tmp * reconVar1_data[512].re));
}

/* End of code generation (FFTImplementationCallback.c) */
