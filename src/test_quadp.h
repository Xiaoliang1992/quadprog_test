/*
 * File: test_quadp.h
 *
 * MATLAB Coder version            : 5.1
 * C/C++ source code generated on  : 02-Nov-2020 23:36:06
 */

#ifndef TEST_QUADP_H
#define TEST_QUADP_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  extern void test_quadp(const double H[9], const double f[3], const double lb[3],
    const double ub[3], double x[3], double *fval);
  extern void test_quadp_initialize(void);
  extern void test_quadp_terminate(void);

#ifdef __cplusplus

}
#endif
#endif

/*
 * File trailer for test_quadp.h
 *
 * [EOF]
 */
