/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * AbstractTimeDomainBeamformer.h
 *
 * Code generation for function 'AbstractTimeDomainBeamformer'
 *
 */

#ifndef ABSTRACTTIMEDOMAINBEAMFORMER_H
#define ABSTRACTTIMEDOMAINBEAMFORMER_H

/* Include files */
#include "only_bf_internal_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void c_AbstractTimeDomainBeamformer_(phased_FrostBeamformer *obj,
                                     const double x[7200],
                                     double x_steered[7200]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (AbstractTimeDomainBeamformer.h) */
