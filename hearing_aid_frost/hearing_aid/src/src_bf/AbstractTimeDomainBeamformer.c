/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * AbstractTimeDomainBeamformer.c
 *
 * Code generation for function 'AbstractTimeDomainBeamformer'
 *
 */

/* Include files */
#include "AbstractTimeDomainBeamformer.h"
#include "FFTImplementationCallback.h"
#include "only_bf_internal_types.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void c_AbstractTimeDomainBeamformer_(phased_FrostBeamformer *obj,
                                     const double x[7200],
                                     double x_steered[7200])
{
  static const double delayFrac[6] = {
      0.49819155242827318,  0.0, -0.49819155242827318,
      -0.49819155242827318, 0.0, 0.49819155242827318};
  static const double dt[6] = {-5.5018084475717268, -0.0, 5.5018084475717268,
                               5.5018084475717268,  -0.0, -5.5018084475717268};
  static double output[7236];
  static const signed char delayInt[6] = {-6, 0, 6, 6, 0, -6};
  int colI;
  int dim;
  int ib;
  int j;
  int k;
  int x_tmp;
  if (obj->cElementDelay.isInitialized != 1) {
    obj->cElementDelay.isSetupComplete = false;
    obj->cElementDelay.isInitialized = 1;
    obj->cElementDelay._pobj1.isInitialized = 0;
    obj->cElementDelay._pobj1.Element = &obj->cElementDelay._pobj0;
    obj->cElementDelay._pobj1.matlabCodegenIsDeleted = false;
    obj->cElementDelay.cSensorArray = &obj->cElementDelay._pobj1;
    obj->cElementDelay.isSetupComplete = true;
  }
  memset(&output[0], 0, 7236U * sizeof(double));
  for (colI = 0; colI < 6; colI++) {
    if (delayFrac[colI] != 0.0) {
      creal_T b_x_data[2048];
      creal_T tmpxd_data[2048];
      double costab_data[1025];
      double sintab_data[1025];
      double costab1q_data[513];
      double sintab_tmp;
      double temp_im;
      double temp_re;
      double temp_re_tmp;
      double twid_re;
      int iDelta2;
      int newStart;
      int orgStart;
      int tmp;
      int xtmp;
      short x_data[2048];
      signed char i;
      costab1q_data[0] = 1.0;
      for (k = 0; k < 256; k++) {
        costab1q_data[k + 1] = cos(0.0030679615757712823 * ((double)k + 1.0));
      }
      for (k = 0; k < 255; k++) {
        costab1q_data[k + 257] =
            sin(0.0030679615757712823 * (512.0 - ((double)k + 257.0)));
      }
      costab1q_data[512] = 0.0;
      costab_data[0] = 1.0;
      sintab_data[0] = 0.0;
      for (k = 0; k < 512; k++) {
        temp_im = costab1q_data[k + 1];
        costab_data[k + 1] = temp_im;
        sintab_tmp = -costab1q_data[511 - k];
        sintab_data[k + 1] = sintab_tmp;
        costab_data[k + 513] = sintab_tmp;
        sintab_data[k + 513] = -temp_im;
      }
      temp_im = dt[colI];
      for (ib = 0; ib < 2048; ib++) {
        x_data[ib] = (short)(ib - 1024);
      }
      for (dim = 0; dim < 2; dim++) {
        orgStart = dim - 1;
        if (dim + 1 <= 1) {
          newStart = 2048;
        } else {
          newStart = 1;
        }
        if (newStart > 1) {
          int vlend2_tmp;
          vlend2_tmp = newStart / 2;
          ib = vlend2_tmp << 1;
          if (ib == newStart) {
            int midoffset;
            int vstride;
            vstride = 1;
            for (k = 0; k <= orgStart; k++) {
              vstride <<= 11;
            }
            midoffset = vlend2_tmp * vstride - 1;
            if (ib == newStart) {
              int i1;
              i1 = 0;
              for (j = 0; j < vstride; j++) {
                i1++;
                ib = i1 + midoffset;
                for (k = 0; k < vlend2_tmp; k++) {
                  orgStart = k * vstride;
                  newStart = (i1 + orgStart) - 1;
                  tmp = x_data[newStart];
                  x_tmp = ib + orgStart;
                  x_data[newStart] = x_data[x_tmp];
                  x_data[x_tmp] = (short)tmp;
                }
              }
            } else {
              int i1;
              i1 = 0;
              for (j = 0; j < vstride; j++) {
                i1++;
                ib = i1 + midoffset;
                xtmp = x_data[ib];
                for (k = 0; k < vlend2_tmp; k++) {
                  orgStart = ib + vstride;
                  x_tmp = (i1 + k * vstride) - 1;
                  x_data[ib] = x_data[x_tmp];
                  x_data[x_tmp] = x_data[orgStart];
                  ib = orgStart;
                }
                x_data[ib] = (short)xtmp;
              }
            }
          } else {
            int i1;
            int midoffset;
            int vstride;
            vstride = 1;
            for (k = 0; k <= orgStart; k++) {
              vstride <<= 11;
            }
            midoffset = vlend2_tmp * vstride - 1;
            i1 = 0;
            orgStart = (newStart - 1) * vstride + -1;
            for (j = 0; j < vstride; j++) {
              i1++;
              newStart = i1 + midoffset;
              ib = (orgStart + j) + 1;
              xtmp = x_data[ib];
              for (k = 0; k < vlend2_tmp; k++) {
                x_tmp = (k + 1) * -vstride;
                iDelta2 = newStart + x_tmp;
                x_data[ib + k * -vstride] = x_data[iDelta2];
                x_data[iDelta2] = x_data[ib + x_tmp];
              }
              x_data[ib + vlend2_tmp * -vstride] = (short)xtmp;
            }
          }
        }
      }
      for (k = 0; k < 2048; k++) {
        sintab_tmp = 6.2831853071795862 * (double)x_data[k] / 2048.0 * -temp_im;
        tmpxd_data[k].re = cos(sintab_tmp);
        tmpxd_data[k].im = sin(sintab_tmp);
        b_x_data[k].re = 0.0;
        b_x_data[k].im = 0.0;
      }
      c_FFTImplementationCallback_doH(&x[1200 * colI], b_x_data, costab_data,
                                      sintab_data);
      for (ib = 0; ib < 2048; ib++) {
        twid_re = b_x_data[ib].re;
        sintab_tmp = tmpxd_data[ib].im;
        temp_im = b_x_data[ib].im;
        temp_re = tmpxd_data[ib].re;
        b_x_data[ib].re = twid_re * temp_re - temp_im * sintab_tmp;
        b_x_data[ib].im = twid_re * sintab_tmp + temp_im * temp_re;
      }
      costab1q_data[0] = 1.0;
      for (k = 0; k < 256; k++) {
        costab1q_data[k + 1] = cos(0.0030679615757712823 * ((double)k + 1.0));
      }
      for (k = 0; k < 255; k++) {
        costab1q_data[k + 257] =
            sin(0.0030679615757712823 * (512.0 - ((double)k + 257.0)));
      }
      costab1q_data[512] = 0.0;
      costab_data[0] = 1.0;
      sintab_data[0] = 0.0;
      for (k = 0; k < 512; k++) {
        temp_im = costab1q_data[k + 1];
        costab_data[k + 1] = temp_im;
        sintab_tmp = costab1q_data[511 - k];
        sintab_data[k + 1] = sintab_tmp;
        costab_data[k + 513] = -sintab_tmp;
        sintab_data[k + 513] = temp_im;
      }
      orgStart = 0;
      newStart = 0;
      for (x_tmp = 0; x_tmp < 2047; x_tmp++) {
        bool tst;
        tmpxd_data[orgStart] = b_x_data[x_tmp];
        orgStart = 2048;
        tst = true;
        while (tst) {
          orgStart >>= 1;
          newStart ^= orgStart;
          tst = ((newStart & orgStart) == 0);
        }
        orgStart = newStart;
      }
      tmpxd_data[orgStart] = b_x_data[2047];
      for (x_tmp = 0; x_tmp <= 2046; x_tmp += 2) {
        temp_re_tmp = tmpxd_data[x_tmp + 1].re;
        sintab_tmp = tmpxd_data[x_tmp + 1].im;
        twid_re = tmpxd_data[x_tmp].re;
        temp_im = tmpxd_data[x_tmp].im;
        tmpxd_data[x_tmp + 1].re = twid_re - temp_re_tmp;
        tmpxd_data[x_tmp + 1].im = temp_im - sintab_tmp;
        tmpxd_data[x_tmp].re = twid_re + temp_re_tmp;
        tmpxd_data[x_tmp].im = temp_im + sintab_tmp;
      }
      ib = 2;
      iDelta2 = 4;
      k = 512;
      xtmp = 2045;
      while (k > 0) {
        for (x_tmp = 0; x_tmp < xtmp; x_tmp += iDelta2) {
          orgStart = x_tmp + ib;
          temp_re = tmpxd_data[orgStart].re;
          temp_im = tmpxd_data[orgStart].im;
          tmpxd_data[orgStart].re = tmpxd_data[x_tmp].re - temp_re;
          tmpxd_data[orgStart].im = tmpxd_data[x_tmp].im - temp_im;
          tmpxd_data[x_tmp].re += temp_re;
          tmpxd_data[x_tmp].im += temp_im;
        }
        orgStart = 1;
        for (j = k; j < 1024; j += k) {
          double twid_im;
          twid_re = costab_data[j];
          twid_im = sintab_data[j];
          x_tmp = orgStart;
          newStart = orgStart + xtmp;
          while (x_tmp < newStart) {
            tmp = x_tmp + ib;
            temp_re_tmp = tmpxd_data[tmp].im;
            sintab_tmp = tmpxd_data[tmp].re;
            temp_re = twid_re * sintab_tmp - twid_im * temp_re_tmp;
            temp_im = twid_re * temp_re_tmp + twid_im * sintab_tmp;
            tmpxd_data[tmp].re = tmpxd_data[x_tmp].re - temp_re;
            tmpxd_data[tmp].im = tmpxd_data[x_tmp].im - temp_im;
            tmpxd_data[x_tmp].re += temp_re;
            tmpxd_data[x_tmp].im += temp_im;
            x_tmp += iDelta2;
          }
          orgStart++;
        }
        k /= 2;
        ib = iDelta2;
        iDelta2 += iDelta2;
        xtmp -= ib;
      }
      for (ib = 0; ib < 2048; ib++) {
        tmpxd_data[ib].re *= 0.00048828125;
        tmpxd_data[ib].im *= 0.00048828125;
      }
      i = delayInt[colI];
      if (i >= 0) {
        orgStart = i;
        newStart = i;
      } else {
        orgStart = 0;
        newStart = 0;
      }
      tmp = i - orgStart;
      for (ib = 0; ib <= tmp + 1199; ib++) {
        output[(newStart + ib) + 1206 * colI] = tmpxd_data[orgStart + ib].re;
      }
    } else {
      int newStart;
      int orgStart;
      int tmp;
      signed char i;
      i = delayInt[colI];
      if (i >= 0) {
        orgStart = 0;
        newStart = i;
      } else {
        orgStart = -i;
        newStart = 0;
      }
      tmp = -orgStart;
      for (ib = 0; ib <= tmp + 1199; ib++) {
        output[(newStart + ib) + 1206 * colI] =
            x[(orgStart + ib) + 1200 * colI];
      }
    }
  }
  for (ib = 0; ib < 6; ib++) {
    memcpy(&x_steered[ib * 1200], &output[ib * 1206], 1200U * sizeof(double));
  }
}

/* End of code generation (AbstractTimeDomainBeamformer.c) */
