/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * BlockSMIWeightsEstimator.h
 *
 * Code generation for function 'BlockSMIWeightsEstimator'
 *
 */

#ifndef BLOCKSMIWEIGHTSESTIMATOR_H
#define BLOCKSMIWEIGHTSESTIMATOR_H

/* Include files */
#include "only_bf_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void c_BlockSMIWeightsEstimator_step(const c_phased_internal_BlockSMIWeigh *obj,
                                     const double x[14400], double w[12]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (BlockSMIWeightsEstimator.h) */
