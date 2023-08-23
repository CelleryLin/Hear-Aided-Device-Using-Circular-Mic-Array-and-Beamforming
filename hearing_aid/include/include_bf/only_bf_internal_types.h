/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * only_bf_internal_types.h
 *
 * Code generation for function 'only_bf'
 *
 */

#ifndef ONLY_BF_INTERNAL_TYPES_H
#define ONLY_BF_INTERNAL_TYPES_H

/* Include files */
#include "only_bf_types.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_c_phased_OmnidirectionalMicroph
#define typedef_c_phased_OmnidirectionalMicroph
typedef struct {
  int __dummy;
} c_phased_OmnidirectionalMicroph;
#endif /* typedef_c_phased_OmnidirectionalMicroph */

#ifndef typedef_phased_UCA
#define typedef_phased_UCA
typedef struct {
  bool matlabCodegenIsDeleted;
  int isInitialized;
  c_phased_OmnidirectionalMicroph *Element;
} phased_UCA;
#endif /* typedef_phased_UCA */

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3
typedef struct {
  unsigned int f1[8];
} cell_wrap_3;
#endif /* typedef_cell_wrap_3 */

#ifndef typedef_phased_ElementDelay
#define typedef_phased_ElementDelay
typedef struct {
  bool matlabCodegenIsDeleted;
  int isInitialized;
  bool isSetupComplete;
  cell_wrap_3 inputVarSize[1];
  phased_UCA *SensorArray;
  phased_UCA *cSensorArray;
  c_phased_OmnidirectionalMicroph _pobj0;
  phased_UCA _pobj1;
} phased_ElementDelay;
#endif /* typedef_phased_ElementDelay */

#ifndef typedef_c_phased_internal_BlockSMIWeigh
#define typedef_c_phased_internal_BlockSMIWeigh
typedef struct {
  bool matlabCodegenIsDeleted;
  int isInitialized;
  bool isSetupComplete;
  bool TunablePropsChanged;
  cell_wrap_3 inputVarSize[1];
  double pNumInputChannels;
  double DiagonalLoadingFactor;
} c_phased_internal_BlockSMIWeigh;
#endif /* typedef_c_phased_internal_BlockSMIWeigh */

#ifndef typedef_phased_FrostBeamformer
#define typedef_phased_FrostBeamformer
typedef struct {
  bool matlabCodegenIsDeleted;
  int isInitialized;
  bool isSetupComplete;
  bool TunablePropsChanged;
  cell_wrap_3 inputVarSize[1];
  double pNumInputChannels;
  phased_UCA *SensorArray;
  phased_ElementDelay cElementDelay;
  double DiagonalLoadingFactor;
  c_phased_internal_BlockSMIWeigh cWeightsEstimator;
  double pDataBuffer[6];
} phased_FrostBeamformer;
#endif /* typedef_phased_FrostBeamformer */

#endif
/* End of code generation (only_bf_internal_types.h) */
