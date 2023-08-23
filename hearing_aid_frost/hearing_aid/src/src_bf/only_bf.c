/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * only_bf.c
 *
 * Code generation for function 'only_bf'
 *
 */

/* Include files */
#include "only_bf.h"
#include "AbstractTimeDomainBeamformer.h"
#include "BlockSMIWeightsEstimator.h"
#include "only_bf_internal_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void only_bf(const double input_sort[7200], double output[1200], phased_UCA uca, phased_FrostBeamformer beamformer)
{
  static double x_presteered[14400];
  static double x_augmented[7206];
  static double x_presteered_temp[7200];
  c_phased_internal_BlockSMIWeigh *obj;
  phased_ElementDelay *b_obj;
  phased_UCA *varargin_2;
  double w[12];
  double val;
  int i;
  int i1;
  int m;
  bool flag;


  /*  BEAMFORMING */
  if (beamformer.isInitialized != 1) {
    beamformer.isInitialized = 1;
    beamformer.pNumInputChannels = 6.0;
    varargin_2 = beamformer.SensorArray;
    beamformer.cElementDelay.isInitialized = 0;
    beamformer.cElementDelay.SensorArray = varargin_2;
    beamformer.cElementDelay.matlabCodegenIsDeleted = false;
    beamformer.cWeightsEstimator.pNumInputChannels = -1.0;
    beamformer.cWeightsEstimator.isInitialized = 0;
    beamformer.cWeightsEstimator.matlabCodegenIsDeleted = false;
    flag = (beamformer.cWeightsEstimator.isInitialized == 1);
    if (flag) {
      beamformer.cWeightsEstimator.TunablePropsChanged = true;
    }
    val = beamformer.DiagonalLoadingFactor;
    beamformer.cWeightsEstimator.DiagonalLoadingFactor = val;
    beamformer.isSetupComplete = true;
    beamformer.TunablePropsChanged = false;
    for (i = 0; i < 6; i++) {
      beamformer.pDataBuffer[i] = 0.0;
    }
  }
  if (beamformer.TunablePropsChanged) {
    beamformer.TunablePropsChanged = false;
    flag = (beamformer.cWeightsEstimator.isInitialized == 1);
    if (flag) {
      beamformer.cWeightsEstimator.TunablePropsChanged = true;
    }
    val = beamformer.DiagonalLoadingFactor;
    beamformer.cWeightsEstimator.DiagonalLoadingFactor = val;
  }
  c_AbstractTimeDomainBeamformer_(&beamformer, input_sort, x_presteered_temp);
  memset(&x_presteered[0], 0, 14400U * sizeof(double));
  for (i = 0; i < 6; i++) {
    x_augmented[1201 * i] = beamformer.pDataBuffer[i];
  }
  for (i = 0; i < 6; i++) {
    memcpy(&x_augmented[i * 1201 + 1], &x_presteered_temp[i * 1200],
           1200U * sizeof(double));
  }
  for (m = 0; m < 2; m++) {
    i = m * 6;
    if (i + 1 > (m + 1) * 6) {
      i = 0;
    }
    for (i1 = 0; i1 < 6; i1++) {
      memcpy(&x_presteered[(i1 + i) * 1200], &x_augmented[(i1 * 1201 - m) + 1],
             1200U * sizeof(double));
    }
  }
  for (i = 0; i < 6; i++) {
    beamformer.pDataBuffer[i] = x_augmented[1201 * i + 1200];
  }
  if (beamformer.cWeightsEstimator.isInitialized != 1) {
    beamformer.cWeightsEstimator.isInitialized = 1;
    beamformer.cWeightsEstimator.pNumInputChannels = 12.0;
    beamformer.cWeightsEstimator.isSetupComplete = true;
    beamformer.cWeightsEstimator.TunablePropsChanged = false;
  }
  if (beamformer.cWeightsEstimator.TunablePropsChanged) {
    beamformer.cWeightsEstimator.TunablePropsChanged = false;
  }
  c_BlockSMIWeightsEstimator_step(&beamformer.cWeightsEstimator, x_presteered,
                                  w);
  for (i = 0; i < 1200; i++) {
    val = 0.0;
    for (i1 = 0; i1 < 12; i1++) {
      val += x_presteered[i + 1200 * i1] * w[i1];
    }
    output[i] = val;
  }
}

/* End of code generation (only_bf.c) */
