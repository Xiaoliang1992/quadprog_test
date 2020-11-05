/*
 * File: test_quadp.c
 *
 * MATLAB Coder version            : 5.1
 * C/C++ source code generated on  : 02-Nov-2020 23:36:06
 */

/* Include Files */
#include "test_quadp.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  double xstar[4];
  double fstar;
  double firstorderopt;
  double lambda[7];
  int state;
  double maxConstr;
  int iterations;
  double searchDir[4];
} struct_T;

#endif                                 /*typedef_struct_T*/

#ifndef typedef_b_struct_T
#define typedef_b_struct_T

typedef struct {
  int MaxIterations;
  double ConstrRelTolFactor;
  double ProbRelTolFactor;
  boolean_T RemainFeasible;
} b_struct_T;

#endif                                 /*typedef_b_struct_T*/

#ifndef typedef_c_struct_T
#define typedef_c_struct_T

typedef struct {
  double workspace_double[28];
  int workspace_int[7];
  int workspace_sort[7];
} c_struct_T;

#endif                                 /*typedef_c_struct_T*/

#ifndef typedef_d_struct_T
#define typedef_d_struct_T

typedef struct {
  int mConstr;
  int mConstrOrig;
  int mConstrMax;
  int nVar;
  int nVarOrig;
  int nVarMax;
  int ldA;
  double lb[4];
  double ub[4];
  int indexLB[4];
  int indexUB[4];
  int indexFixed[4];
  int mEqRemoved;
  double ATwset[28];
  double bwset[7];
  int nActiveConstr;
  double maxConstrWorkspace[7];
  int sizes[5];
  int sizesNormal[5];
  int sizesPhaseOne[5];
  int sizesRegularized[5];
  int sizesRegPhaseOne[5];
  int isActiveIdx[6];
  int isActiveIdxNormal[6];
  int isActiveIdxPhaseOne[6];
  int isActiveIdxRegularized[6];
  int isActiveIdxRegPhaseOne[6];
  boolean_T isActiveConstr[7];
  int Wid[7];
  int Wlocalidx[7];
  int nWConstr[5];
  int probType;
  double SLACK0;
} d_struct_T;

#endif                                 /*typedef_d_struct_T*/

#ifndef typedef_e_struct_T
#define typedef_e_struct_T

typedef struct {
  int ldq;
  double QR[49];
  double Q[49];
  int jpvt[7];
  int mrows;
  int ncols;
  double tau[7];
  int minRowCol;
  boolean_T usedPivoting;
} e_struct_T;

#endif                                 /*typedef_e_struct_T*/

#ifndef typedef_f_struct_T
#define typedef_f_struct_T

typedef struct {
  double FMat[49];
  int ldm;
  int ndims;
  int info;
  double scaleFactor;
  boolean_T ConvexCheck;
  double regTol_;
  double workspace_[336];
  double workspace2_[336];
} f_struct_T;

#endif                                 /*typedef_f_struct_T*/

#ifndef typedef_g_struct_T
#define typedef_g_struct_T

typedef struct {
  double grad[4];
  double Hx[3];
  boolean_T hasLinear;
  int nvar;
  int maxVar;
  double beta;
  double rho;
  int objtype;
  int prev_objtype;
  int prev_nvar;
  boolean_T prev_hasLinear;
  double gammaScalar;
} g_struct_T;

#endif                                 /*typedef_g_struct_T*/

#ifndef typedef_h_struct_T
#define typedef_h_struct_T

typedef struct {
  double InitDamping;
  char FiniteDifferenceType[7];
  boolean_T SpecifyObjectiveGradient;
  boolean_T ScaleProblem;
  boolean_T SpecifyConstraintGradient;
  boolean_T NonFiniteSupport;
  boolean_T IterDisplaySQP;
  double FiniteDifferenceStepSize;
  double MaxFunctionEvaluations;
  boolean_T IterDisplayQP;
  double PricingTolerance;
  char Algorithm[10];
  double ObjectiveLimit;
  double ConstraintTolerance;
  double OptimalityTolerance;
  double StepTolerance;
  double MaxIterations;
  double FunctionTolerance;
  char SolverName[8];
  boolean_T CheckGradients;
  char Diagnostics[3];
  double DiffMaxChange;
  double DiffMinChange;
  char Display[5];
  char FunValCheck[3];
  boolean_T UseParallel;
  char LinearSolver[4];
  char SubproblemAlgorithm[2];
} h_struct_T;

#endif                                 /*typedef_h_struct_T*/

/* Variable Definitions */
static boolean_T isInitialized_test_quadp = false;

/* Function Declarations */
static void PresolveWorkingSet(struct_T *solution, c_struct_T *memspace,
  d_struct_T *workingset, e_struct_T *qrmanager);
static void RemoveDependentIneq_(d_struct_T *workingset, e_struct_T *qrmanager,
  c_struct_T *memspace, double tolfactor);
static void addBoundToActiveSetMatrix_(d_struct_T *obj, int TYPE, int idx_local);
static double b_maxConstraintViolation(const d_struct_T *obj, const double x[4]);
static void b_xgemm(int m, int n, int k, const double A[49], int ia0, const
                    double B[28], double C[49]);
static double b_xnrm2(int n, const double x[4]);
static void checkStoppingAndUpdateFval(int *activeSetChangeID, const double f[3],
  struct_T *solution, c_struct_T *memspace, const g_struct_T *objective, const
  d_struct_T *workingset, e_struct_T *qrmanager, double options_ObjectiveLimit,
  int runTimeOptions_MaxIterations, boolean_T updateFval);
static void computeFirstOrderOpt(struct_T *solution, const g_struct_T *objective,
  int workingset_nVar, const double workingset_ATwset[28], int
  workingset_nActiveConstr, double workspace[28]);
static double computeFval(const g_struct_T *obj, double workspace[28], const
  double H[9], const double f[3], const double x[4]);
static double computeFval_ReuseHx(const g_struct_T *obj, double workspace[28],
  const double f[3], const double x[4]);
static void computeGrad_StoreHx(g_struct_T *obj, const double H[9], const double
  f[3], const double x[4]);
static void computeQ_(e_struct_T *obj, int nrows);
static void compute_deltax(const double H[9], struct_T *solution, c_struct_T
  *memspace, const e_struct_T *qrmanager, f_struct_T *cholmanager, const
  g_struct_T *objective);
static void compute_lambda(double workspace[28], struct_T *solution, const
  g_struct_T *objective, const e_struct_T *qrmanager);
static void countsort(int x[7], int xLen, int workspace[7], int xMin, int xMax);
static void deleteColMoveEnd(e_struct_T *obj, int idx);
static void driver(const double H[9], const double f[3], struct_T *solution,
                   c_struct_T *memspace, d_struct_T *workingset, b_struct_T
                   runTimeOptions, e_struct_T *qrmanager, f_struct_T
                   *cholmanager, g_struct_T *objective);
static void factor(f_struct_T *obj, const double A[9], int ndims, int ldA);
static void factorQR(e_struct_T *obj, const double A[28], int mrows, int ncols);
static void factorQRE(e_struct_T *obj, const double A[28], int mrows, int ncols);
static boolean_T feasibleX0ForWorkingSet(double workspace[28], double xCurrent[4],
  const d_struct_T *workingset, e_struct_T *qrmanager);
static void feasibleratiotest(const double solution_xstar[4], const double
  solution_searchDir[4], int workingset_nVar, const double workingset_lb[4],
  const double workingset_ub[4], const int workingset_indexLB[4], const int
  workingset_indexUB[4], const int workingset_sizes[5], const int
  workingset_isActiveIdx[6], const boolean_T workingset_isActiveConstr[7], const
  int workingset_nWConstr[5], boolean_T isPhaseOne, double *alpha, boolean_T
  *newBlocking, int *constrType, int *constrIdx);
static void fullColLDL2_(f_struct_T *obj, int LD_offset, int NColsRemain, double
  REG_PRIMAL);
static void iterate(const double H[9], const double f[3], struct_T *solution,
                    c_struct_T *memspace, d_struct_T *workingset, e_struct_T
                    *qrmanager, f_struct_T *cholmanager, g_struct_T *objective,
                    double options_ObjectiveLimit, double options_StepTolerance,
                    int runTimeOptions_MaxIterations, double
                    runTimeOptions_ProbRelTolFactor, boolean_T
                    runTimeOptions_RemainFeasible);
static void linearForm_(boolean_T obj_hasLinear, int obj_nvar, double workspace
  [28], const double H[9], const double f[3], const double x[4]);
static double maxConstraintViolation(const d_struct_T *obj, const double x[28],
  int ix0);
static void partialColLDL3_(f_struct_T *obj, int LD_offset, int NColsRemain,
  double REG_PRIMAL);
static void phaseone(const double H[9], const double f[3], struct_T *solution,
                     c_struct_T *memspace, d_struct_T *workingset, e_struct_T
                     *qrmanager, h_struct_T *options, const b_struct_T
                     *runTimeOptions, f_struct_T *cholmanager, g_struct_T
                     *objective);
static void qrf(double A[49], int m, int n, int nfxd, double tau[7]);
static void removeConstr(d_struct_T *obj, int idx_global);
static void removeEqConstr(d_struct_T *obj, int idx_global);
static double rt_hypotd_snf(double u0, double u1);
static void setProblemType(d_struct_T *obj, int PROBLEM_TYPE);
static void solve(const f_struct_T *obj, double rhs[4]);
static void squareQ_appendCol(e_struct_T *obj, const double vec[28], int iv0);
static void xgemm(int m, int n, int k, const double A[9], int lda, const double
                  B[49], int ib0, double C[28]);
static void xgeqp3(double A[49], int m, int n, int jpvt[7], double tau[7]);
static double xnrm2(int n, const double x[49], int ix0);
static void xrotg(double *a, double *b, double *c, double *s);
static void xzlarf(int m, int n, int iv0, double tau, double C[49], int ic0,
                   double work[7]);
static double xzlarfg(int n, double *alpha1, double x[49], int ix0);

/* Function Definitions */
/*
 * Arguments    : struct_T *solution
 *                c_struct_T *memspace
 *                d_struct_T *workingset
 *                e_struct_T *qrmanager
 * Return Type  : void
 */
static void PresolveWorkingSet(struct_T *solution, c_struct_T *memspace,
  d_struct_T *workingset, e_struct_T *qrmanager)
{
  double qtb;
  double tol;
  int idx;
  int idx_col;
  int ix;
  int k;
  int mTotalWorkingEq_tmp_tmp;
  int mWorkingFixed;
  int nVar;
  boolean_T exitg1;
  boolean_T guard1 = false;
  boolean_T okWorkingSet;
  solution->state = 82;
  nVar = workingset->nVar - 1;
  mWorkingFixed = workingset->nWConstr[0];
  mTotalWorkingEq_tmp_tmp = workingset->nWConstr[1] + workingset->nWConstr[0];
  idx_col = 0;
  if (mTotalWorkingEq_tmp_tmp > 0) {
    for (ix = 0; ix < mTotalWorkingEq_tmp_tmp; ix++) {
      for (idx_col = 0; idx_col <= nVar; idx_col++) {
        qrmanager->QR[ix + 7 * idx_col] = workingset->ATwset[idx_col + (ix << 2)];
      }
    }

    idx_col = mTotalWorkingEq_tmp_tmp - workingset->nVar;
    if (0 > idx_col) {
      idx_col = 0;
    }

    if (0 <= nVar) {
      memset(&qrmanager->jpvt[0], 0, (nVar + 1) * sizeof(int));
    }

    qrmanager->usedPivoting = true;
    qrmanager->mrows = mTotalWorkingEq_tmp_tmp;
    qrmanager->ncols = workingset->nVar;
    nVar = workingset->nVar;
    if (mTotalWorkingEq_tmp_tmp < nVar) {
      nVar = mTotalWorkingEq_tmp_tmp;
    }

    qrmanager->minRowCol = nVar;
    xgeqp3(qrmanager->QR, mTotalWorkingEq_tmp_tmp, workingset->nVar,
           qrmanager->jpvt, qrmanager->tau);
    tol = 100.0 * (double)workingset->nVar * 2.2204460492503131E-16;
    nVar = workingset->nVar;
    if (nVar >= mTotalWorkingEq_tmp_tmp) {
      nVar = mTotalWorkingEq_tmp_tmp;
    }

    while ((nVar > 0) && (fabs(qrmanager->QR[(nVar + 7 * (nVar - 1)) - 1]) < tol))
    {
      nVar--;
      idx_col++;
    }

    if (idx_col > 0) {
      computeQ_(qrmanager, mTotalWorkingEq_tmp_tmp);
      idx = 0;
      exitg1 = false;
      while ((!exitg1) && (idx <= idx_col - 1)) {
        ix = 7 * ((mTotalWorkingEq_tmp_tmp - idx) - 1);
        qtb = 0.0;
        nVar = 0;
        for (k = 0; k < mTotalWorkingEq_tmp_tmp; k++) {
          qtb += qrmanager->Q[ix] * workingset->bwset[nVar];
          ix++;
          nVar++;
        }

        if (fabs(qtb) >= tol) {
          idx_col = -1;
          exitg1 = true;
        } else {
          idx++;
        }
      }
    }

    if (idx_col > 0) {
      for (idx = 0; idx < mWorkingFixed; idx++) {
        qrmanager->jpvt[idx] = 1;
      }

      nVar = workingset->nWConstr[0] + 1;
      if (nVar <= mTotalWorkingEq_tmp_tmp) {
        memset(&qrmanager->jpvt[nVar + -1], 0, ((mTotalWorkingEq_tmp_tmp - nVar)
                + 1) * sizeof(int));
      }

      factorQRE(qrmanager, workingset->ATwset, workingset->nVar,
                mTotalWorkingEq_tmp_tmp);
      for (idx = 0; idx < idx_col; idx++) {
        memspace->workspace_int[idx] = qrmanager->jpvt[(mTotalWorkingEq_tmp_tmp
          - idx_col) + idx];
      }

      countsort(memspace->workspace_int, idx_col, memspace->workspace_sort, 1,
                mTotalWorkingEq_tmp_tmp);
      for (idx = idx_col; idx >= 1; idx--) {
        removeEqConstr(workingset, memspace->workspace_int[idx - 1]);
      }
    }
  }

  if (idx_col != -1) {
    RemoveDependentIneq_(workingset, qrmanager, memspace, 100.0);
    okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_double,
      solution->xstar, workingset, qrmanager);
    guard1 = false;
    if (!okWorkingSet) {
      RemoveDependentIneq_(workingset, qrmanager, memspace, 1000.0);
      okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_double,
        solution->xstar, workingset, qrmanager);
      if (!okWorkingSet) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 && (workingset->nWConstr[0] + workingset->nWConstr[1] ==
                   workingset->nVar)) {
      tol = b_maxConstraintViolation(workingset, solution->xstar);
      if (tol > 1.0E-8) {
        solution->state = -2;
      }
    }
  } else {
    solution->state = -3;
    nVar = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
    ix = workingset->nActiveConstr;
    for (idx_col = nVar; idx_col <= ix; idx_col++) {
      workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
        Wid[idx_col - 1] - 1] + workingset->Wlocalidx[idx_col - 1]) - 2] = false;
    }

    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  }
}

/*
 * Arguments    : d_struct_T *workingset
 *                e_struct_T *qrmanager
 *                c_struct_T *memspace
 *                double tolfactor
 * Return Type  : void
 */
static void RemoveDependentIneq_(d_struct_T *workingset, e_struct_T *qrmanager,
  c_struct_T *memspace, double tolfactor)
{
  double tol;
  int b_idx;
  int i;
  int idx;
  int nActiveConstr;
  int nDepIneq;
  int nFixedConstr;
  nActiveConstr = workingset->nActiveConstr;
  nFixedConstr = workingset->nWConstr[1] + workingset->nWConstr[0];
  if ((workingset->nWConstr[2] + workingset->nWConstr[3]) + workingset->
      nWConstr[4] > 0) {
    tol = tolfactor * (double)workingset->nVar * 2.2204460492503131E-16;
    for (idx = 0; idx < nFixedConstr; idx++) {
      qrmanager->jpvt[idx] = 1;
    }

    i = nFixedConstr + 1;
    if (i <= nActiveConstr) {
      memset(&qrmanager->jpvt[i + -1], 0, ((nActiveConstr - i) + 1) * sizeof(int));
    }

    factorQRE(qrmanager, workingset->ATwset, workingset->nVar,
              workingset->nActiveConstr);
    nDepIneq = 0;
    for (idx = workingset->nActiveConstr; idx > workingset->nVar; idx--) {
      nDepIneq++;
      memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx - 1];
    }

    if (idx <= workingset->nVar) {
      while ((idx > nFixedConstr) && (fabs(qrmanager->QR[(idx + 7 * (idx - 1)) -
               1]) < tol)) {
        nDepIneq++;
        memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx - 1];
        idx--;
      }
    }

    countsort(memspace->workspace_int, nDepIneq, memspace->workspace_sort,
              nFixedConstr + 1, workingset->nActiveConstr);
    for (idx = nDepIneq; idx >= 1; idx--) {
      nActiveConstr = memspace->workspace_int[idx - 1] - 1;
      nFixedConstr = workingset->Wid[nActiveConstr] - 1;
      workingset->isActiveConstr[(workingset->isActiveIdx[nFixedConstr] +
        workingset->Wlocalidx[nActiveConstr]) - 2] = false;
      workingset->Wid[nActiveConstr] = workingset->Wid[workingset->nActiveConstr
        - 1];
      workingset->Wlocalidx[nActiveConstr] = workingset->Wlocalidx
        [workingset->nActiveConstr - 1];
      i = workingset->nVar;
      for (b_idx = 0; b_idx < i; b_idx++) {
        workingset->ATwset[b_idx + (nActiveConstr << 2)] = workingset->
          ATwset[b_idx + ((workingset->nActiveConstr - 1) << 2)];
      }

      workingset->bwset[nActiveConstr] = workingset->bwset
        [workingset->nActiveConstr - 1];
      workingset->nActiveConstr--;
      workingset->nWConstr[nFixedConstr]--;
    }
  }
}

/*
 * Arguments    : d_struct_T *obj
 *                int TYPE
 *                int idx_local
 * Return Type  : void
 */
static void addBoundToActiveSetMatrix_(d_struct_T *obj, int TYPE, int idx_local)
{
  int i;
  int i1;
  int idx_bnd_local;
  obj->nWConstr[TYPE - 1]++;
  obj->isActiveConstr[(obj->isActiveIdx[TYPE - 1] + idx_local) - 2] = true;
  obj->nActiveConstr++;
  obj->Wid[obj->nActiveConstr - 1] = TYPE;
  obj->Wlocalidx[obj->nActiveConstr - 1] = idx_local;
  i = obj->nActiveConstr - 1;
  switch (TYPE) {
   case 5:
    idx_bnd_local = obj->indexUB[idx_local - 1];
    obj->bwset[obj->nActiveConstr - 1] = obj->ub[idx_bnd_local - 1];
    break;

   default:
    idx_bnd_local = obj->indexLB[idx_local - 1];
    obj->bwset[obj->nActiveConstr - 1] = obj->lb[idx_bnd_local - 1];
    break;
  }

  if (0 <= idx_bnd_local - 2) {
    memset(&obj->ATwset[i * 4], 0, (((idx_bnd_local + i * 4) - i * 4) + -1) *
           sizeof(double));
  }

  obj->ATwset[(idx_bnd_local + ((obj->nActiveConstr - 1) << 2)) - 1] = 2.0 *
    (double)(TYPE == 5) - 1.0;
  idx_bnd_local++;
  i1 = obj->nVar;
  if (idx_bnd_local <= i1) {
    memset(&obj->ATwset[(idx_bnd_local + i * 4) + -1], 0, ((((i1 + i * 4) -
              idx_bnd_local) - i * 4) + 1) * sizeof(double));
  }

  switch (obj->probType) {
   case 3:
   case 2:
    break;

   default:
    obj->ATwset[(obj->nVar + ((obj->nActiveConstr - 1) << 2)) - 1] = -1.0;
    break;
  }
}

/*
 * Arguments    : const d_struct_T *obj
 *                const double x[4]
 * Return Type  : double
 */
static double b_maxConstraintViolation(const d_struct_T *obj, const double x[4])
{
  double u1;
  double v;
  int idx;
  int mFixed;
  int mLB;
  int mUB;
  mLB = obj->sizes[3];
  mUB = obj->sizes[4];
  mFixed = obj->sizes[0];
  v = 0.0;
  if (obj->sizes[3] > 0) {
    for (idx = 0; idx < mLB; idx++) {
      u1 = -x[obj->indexLB[idx] - 1] - obj->lb[obj->indexLB[idx] - 1];
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    for (idx = 0; idx < mUB; idx++) {
      u1 = x[obj->indexUB[idx] - 1] - obj->ub[obj->indexUB[idx] - 1];
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    for (idx = 0; idx < mFixed; idx++) {
      mLB = obj->indexFixed[idx] - 1;
      u1 = fabs(x[mLB] - obj->ub[mLB]);
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  return v;
}

/*
 * Arguments    : int m
 *                int n
 *                int k
 *                const double A[49]
 *                int ia0
 *                const double B[28]
 *                double C[49]
 * Return Type  : void
 */
static void b_xgemm(int m, int n, int k, const double A[49], int ia0, const
                    double B[28], double C[49])
{
  double temp;
  int ar;
  int br;
  int cr;
  int i;
  int i1;
  int ic;
  int lastColC;
  int w;
  if ((m != 0) && (n != 0)) {
    lastColC = 7 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 7) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
      }
    }

    br = -1;
    for (cr = 0; cr <= lastColC; cr += 7) {
      ar = ia0;
      i = cr + 1;
      i1 = cr + m;
      for (ic = i; ic <= i1; ic++) {
        temp = 0.0;
        for (w = 0; w < k; w++) {
          temp += A[(w + ar) - 1] * B[(w + br) + 1];
        }

        C[ic - 1] += temp;
        ar += 7;
      }

      br += 7;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x[4]
 * Return Type  : double
 */
static double b_xnrm2(int n, const double x[4])
{
  double absxk;
  double scale;
  double t;
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[0]);
    } else {
      scale = 3.3121686421112381E-170;
      for (k = 0; k < n; k++) {
        absxk = fabs(x[k]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : int *activeSetChangeID
 *                const double f[3]
 *                struct_T *solution
 *                c_struct_T *memspace
 *                const g_struct_T *objective
 *                const d_struct_T *workingset
 *                e_struct_T *qrmanager
 *                double options_ObjectiveLimit
 *                int runTimeOptions_MaxIterations
 *                boolean_T updateFval
 * Return Type  : void
 */
static void checkStoppingAndUpdateFval(int *activeSetChangeID, const double f[3],
  struct_T *solution, c_struct_T *memspace, const g_struct_T *objective, const
  d_struct_T *workingset, e_struct_T *qrmanager, double options_ObjectiveLimit,
  int runTimeOptions_MaxIterations, boolean_T updateFval)
{
  double constrViolation_new;
  int nVar;
  boolean_T nonDegenerateWset;
  solution->iterations++;
  nVar = objective->nvar - 1;
  if ((solution->iterations >= runTimeOptions_MaxIterations) &&
      ((solution->state != 1) || (objective->objtype == 5))) {
    solution->state = 0;
  }

  if (solution->iterations - solution->iterations / 50 * 50 == 0) {
    solution->maxConstr = b_maxConstraintViolation(workingset, solution->xstar);
    if (solution->maxConstr > 1.0E-8) {
      if (0 <= nVar) {
        memcpy(&solution->searchDir[0], &solution->xstar[0], (nVar + 1) * sizeof
               (double));
      }

      nonDegenerateWset = feasibleX0ForWorkingSet(memspace->workspace_double,
        solution->searchDir, workingset, qrmanager);
      if ((!nonDegenerateWset) && (solution->state != 0)) {
        solution->state = -2;
      }

      *activeSetChangeID = 0;
      constrViolation_new = b_maxConstraintViolation(workingset,
        solution->searchDir);
      if (constrViolation_new < solution->maxConstr) {
        if (0 <= nVar) {
          memcpy(&solution->xstar[0], &solution->searchDir[0], (nVar + 1) *
                 sizeof(double));
        }

        solution->maxConstr = constrViolation_new;
      }
    }
  }

  if (updateFval) {
    solution->fstar = computeFval_ReuseHx(objective, memspace->workspace_double,
      f, solution->xstar);
    if ((solution->fstar < options_ObjectiveLimit) && ((solution->state != 0) ||
         (objective->objtype != 5))) {
      solution->state = 2;
    }
  }
}

/*
 * Arguments    : struct_T *solution
 *                const g_struct_T *objective
 *                int workingset_nVar
 *                const double workingset_ATwset[28]
 *                int workingset_nActiveConstr
 *                double workspace[28]
 * Return Type  : void
 */
static void computeFirstOrderOpt(struct_T *solution, const g_struct_T *objective,
  int workingset_nVar, const double workingset_ATwset[28], int
  workingset_nActiveConstr, double workspace[28])
{
  double s;
  double smax;
  int i;
  int ia;
  int iac;
  int idxmax;
  int ix;
  int iy;
  if (0 <= workingset_nVar - 1) {
    memcpy(&workspace[0], &objective->grad[0], workingset_nVar * sizeof(double));
  }

  if (workingset_nActiveConstr != 0) {
    ix = 0;
    idxmax = ((workingset_nActiveConstr - 1) << 2) + 1;
    for (iac = 1; iac <= idxmax; iac += 4) {
      iy = 0;
      i = (iac + workingset_nVar) - 1;
      for (ia = iac; ia <= i; ia++) {
        workspace[iy] += workingset_ATwset[ia - 1] * solution->lambda[ix];
        iy++;
      }

      ix++;
    }
  }

  idxmax = 1;
  ix = 0;
  smax = fabs(workspace[0]);
  for (iy = 2; iy <= workingset_nVar; iy++) {
    ix++;
    s = fabs(workspace[ix]);
    if (s > smax) {
      idxmax = iy;
      smax = s;
    }
  }

  solution->firstorderopt = fabs(workspace[idxmax - 1]);
}

/*
 * Arguments    : const g_struct_T *obj
 *                double workspace[28]
 *                const double H[9]
 *                const double f[3]
 *                const double x[4]
 * Return Type  : double
 */
static double computeFval(const g_struct_T *obj, double workspace[28], const
  double H[9], const double f[3], const double x[4])
{
  double val;
  int idx;
  int ixlast;
  switch (obj->objtype) {
   case 5:
    val = obj->gammaScalar * x[obj->nvar - 1];
    break;

   case 3:
    linearForm_(obj->hasLinear, obj->nvar, workspace, H, f, x);
    val = 0.0;
    if (obj->nvar >= 1) {
      ixlast = obj->nvar;
      for (idx = 0; idx < ixlast; idx++) {
        val += x[idx] * workspace[idx];
      }
    }
    break;

   default:
    linearForm_(obj->hasLinear, obj->nvar, workspace, H, f, x);
    ixlast = obj->nvar + 1;
    for (idx = ixlast; idx < 4; idx++) {
      workspace[idx - 1] = 0.0 * x[idx - 1];
    }

    val = (x[0] * workspace[0] + x[1] * workspace[1]) + x[2] * workspace[2];
    break;
  }

  return val;
}

/*
 * Arguments    : const g_struct_T *obj
 *                double workspace[28]
 *                const double f[3]
 *                const double x[4]
 * Return Type  : double
 */
static double computeFval_ReuseHx(const g_struct_T *obj, double workspace[28],
  const double f[3], const double x[4])
{
  double val;
  int idx;
  int ixlast;
  switch (obj->objtype) {
   case 5:
    val = obj->gammaScalar * x[obj->nvar - 1];
    break;

   case 3:
    if (obj->hasLinear) {
      ixlast = obj->nvar;
      for (idx = 0; idx < ixlast; idx++) {
        workspace[idx] = 0.5 * obj->Hx[idx] + f[idx];
      }

      val = 0.0;
      if (obj->nvar >= 1) {
        ixlast = obj->nvar;
        for (idx = 0; idx < ixlast; idx++) {
          val += x[idx] * workspace[idx];
        }
      }
    } else {
      val = 0.0;
      if (obj->nvar >= 1) {
        ixlast = obj->nvar;
        for (idx = 0; idx < ixlast; idx++) {
          val += x[idx] * obj->Hx[idx];
        }
      }

      val *= 0.5;
    }
    break;

   default:
    if (obj->hasLinear) {
      ixlast = obj->nvar;
      if (0 <= ixlast - 1) {
        memcpy(&workspace[0], &f[0], ixlast * sizeof(double));
      }

      ixlast = 2 - obj->nvar;
      for (idx = 0; idx <= ixlast; idx++) {
        workspace[obj->nvar + idx] = 0.0;
      }

      workspace[0] += 0.5 * obj->Hx[0];
      workspace[1] += 0.5 * obj->Hx[1];
      workspace[2] += 0.5 * obj->Hx[2];
      val = (x[0] * workspace[0] + x[1] * workspace[1]) + x[2] * workspace[2];
    } else {
      val = 0.5 * ((x[0] * obj->Hx[0] + x[1] * obj->Hx[1]) + x[2] * obj->Hx[2]);
      ixlast = obj->nvar + 1;
      for (idx = ixlast; idx < 4; idx++) {
        val += x[idx - 1] * 0.0;
      }
    }
    break;
  }

  return val;
}

/*
 * Arguments    : g_struct_T *obj
 *                const double H[9]
 *                const double f[3]
 *                const double x[4]
 * Return Type  : void
 */
static void computeGrad_StoreHx(g_struct_T *obj, const double H[9], const double
  f[3], const double x[4])
{
  int i;
  int i1;
  int ia;
  int iac;
  int ix;
  int ixlast;
  int iy;
  int lda;
  switch (obj->objtype) {
   case 5:
    i = obj->nvar;
    if (0 <= i - 2) {
      memset(&obj->grad[0], 0, (i + -1) * sizeof(double));
    }

    obj->grad[obj->nvar - 1] = obj->gammaScalar;
    break;

   case 3:
    ixlast = obj->nvar - 1;
    lda = obj->nvar;
    if (obj->nvar != 0) {
      if (0 <= ixlast) {
        memset(&obj->Hx[0], 0, (ixlast + 1) * sizeof(double));
      }

      ix = 0;
      i = obj->nvar * (obj->nvar - 1) + 1;
      for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
        iy = 0;
        i1 = iac + ixlast;
        for (ia = iac; ia <= i1; ia++) {
          obj->Hx[iy] += H[ia - 1] * x[ix];
          iy++;
        }

        ix++;
      }
    }

    i = obj->nvar;
    if (0 <= i - 1) {
      memcpy(&obj->grad[0], &obj->Hx[0], i * sizeof(double));
    }

    if (obj->hasLinear && (obj->nvar >= 1)) {
      ixlast = obj->nvar - 1;
      for (lda = 0; lda <= ixlast; lda++) {
        obj->grad[lda] += f[lda];
      }
    }
    break;

   default:
    ixlast = obj->nvar - 1;
    lda = obj->nvar;
    if (obj->nvar != 0) {
      if (0 <= ixlast) {
        memset(&obj->Hx[0], 0, (ixlast + 1) * sizeof(double));
      }

      ix = 0;
      i = obj->nvar * (obj->nvar - 1) + 1;
      for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
        iy = 0;
        i1 = iac + ixlast;
        for (ia = iac; ia <= i1; ia++) {
          obj->Hx[iy] += H[ia - 1] * x[ix];
          iy++;
        }

        ix++;
      }
    }

    i = obj->nvar + 1;
    for (ixlast = i; ixlast < 4; ixlast++) {
      obj->Hx[ixlast - 1] = 0.0 * x[ixlast - 1];
    }

    obj->grad[0] = obj->Hx[0];
    obj->grad[1] = obj->Hx[1];
    obj->grad[2] = obj->Hx[2];
    if (obj->hasLinear && (obj->nvar >= 1)) {
      ixlast = obj->nvar - 1;
      for (lda = 0; lda <= ixlast; lda++) {
        obj->grad[lda] += f[lda];
      }
    }
    break;
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                int nrows
 * Return Type  : void
 */
static void computeQ_(e_struct_T *obj, int nrows)
{
  double work[7];
  double c;
  int b_i;
  int exitg1;
  int i;
  int i1;
  int iQR0;
  int ia;
  int iaii;
  int idx;
  int itau;
  int ix;
  int j;
  int jy;
  int lastc;
  int lastv;
  int m;
  boolean_T exitg2;
  i = obj->minRowCol;
  for (idx = 0; idx < i; idx++) {
    iQR0 = 7 * idx + idx;
    jy = obj->mrows - idx;
    if (0 <= jy - 2) {
      memcpy(&obj->Q[iQR0 + 1], &obj->QR[iQR0 + 1], (((jy + iQR0) - iQR0) + -1) *
             sizeof(double));
    }
  }

  m = obj->mrows;
  jy = obj->minRowCol;
  if (nrows >= 1) {
    i = nrows - 1;
    for (j = jy; j <= i; j++) {
      ia = j * 7;
      i1 = m - 1;
      if (0 <= i1) {
        memset(&obj->Q[ia], 0, (((i1 + ia) - ia) + 1) * sizeof(double));
      }

      obj->Q[ia + j] = 1.0;
    }

    itau = obj->minRowCol - 1;
    for (b_i = 0; b_i < 7; b_i++) {
      work[b_i] = 0.0;
    }

    for (b_i = obj->minRowCol; b_i >= 1; b_i--) {
      iaii = b_i + (b_i - 1) * 7;
      if (b_i < nrows) {
        obj->Q[iaii - 1] = 1.0;
        jy = iaii + 7;
        if (obj->tau[itau] != 0.0) {
          lastv = m - b_i;
          iQR0 = (iaii + m) - b_i;
          while ((lastv + 1 > 0) && (obj->Q[iQR0 - 1] == 0.0)) {
            lastv--;
            iQR0--;
          }

          lastc = (nrows - b_i) - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc + 1 > 0)) {
            iQR0 = (iaii + lastc * 7) + 7;
            ia = iQR0;
            do {
              exitg1 = 0;
              if (ia <= iQR0 + lastv) {
                if (obj->Q[ia - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  ia++;
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
        } else {
          lastv = -1;
          lastc = -1;
        }

        if (lastv + 1 > 0) {
          if (lastc + 1 != 0) {
            if (0 <= lastc) {
              memset(&work[0], 0, (lastc + 1) * sizeof(double));
            }

            iQR0 = 0;
            i = (iaii + 7 * lastc) + 7;
            for (idx = jy; idx <= i; idx += 7) {
              ix = iaii;
              c = 0.0;
              i1 = idx + lastv;
              for (ia = idx; ia <= i1; ia++) {
                c += obj->Q[ia - 1] * obj->Q[ix - 1];
                ix++;
              }

              work[iQR0] += c;
              iQR0++;
            }
          }

          if (!(-obj->tau[itau] == 0.0)) {
            iQR0 = iaii;
            jy = 0;
            for (j = 0; j <= lastc; j++) {
              if (work[jy] != 0.0) {
                c = work[jy] * -obj->tau[itau];
                ix = iaii;
                i = iQR0 + 7;
                i1 = lastv + iQR0;
                for (idx = i; idx <= i1 + 7; idx++) {
                  obj->Q[idx - 1] += obj->Q[ix - 1] * c;
                  ix++;
                }
              }

              jy++;
              iQR0 += 7;
            }
          }
        }
      }

      if (b_i < m) {
        iQR0 = iaii + 1;
        i = (iaii + m) - b_i;
        for (jy = iQR0; jy <= i; jy++) {
          obj->Q[jy - 1] *= -obj->tau[itau];
        }
      }

      obj->Q[iaii - 1] = 1.0 - obj->tau[itau];
      for (j = 0; j <= b_i - 2; j++) {
        obj->Q[(iaii - j) - 2] = 0.0;
      }

      itau--;
    }
  }
}

/*
 * Arguments    : const double H[9]
 *                struct_T *solution
 *                c_struct_T *memspace
 *                const e_struct_T *qrmanager
 *                f_struct_T *cholmanager
 *                const g_struct_T *objective
 * Return Type  : void
 */
static void compute_deltax(const double H[9], struct_T *solution, c_struct_T
  *memspace, const e_struct_T *qrmanager, f_struct_T *cholmanager, const
  g_struct_T *objective)
{
  double SCALED_REG_PRIMAL;
  double s;
  double smax;
  int exitg2;
  int i;
  int ia;
  int iac;
  int idx_row;
  int ix;
  int mNull_tmp;
  int nVar;
  int nVars;
  int nullStartIdx_tmp;
  int order;
  boolean_T exitg1;
  nVar = qrmanager->mrows - 1;
  mNull_tmp = qrmanager->mrows - qrmanager->ncols;
  if (mNull_tmp <= 0) {
    if (0 <= nVar) {
      memset(&solution->searchDir[0], 0, (nVar + 1) * sizeof(double));
    }
  } else {
    for (order = 0; order <= nVar; order++) {
      solution->searchDir[order] = -objective->grad[order];
    }

    if (qrmanager->ncols <= 0) {
      switch (objective->objtype) {
       case 5:
        break;

       case 3:
        factor(cholmanager, H, qrmanager->mrows, qrmanager->mrows);
        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          solve(cholmanager, solution->searchDir);
        }
        break;

       default:
        factor(cholmanager, H, objective->nvar, objective->nvar);
        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          solve(cholmanager, solution->searchDir);
          nVars = objective->nvar + 1;
          i = qrmanager->mrows;
          for (idx_row = nVars; idx_row <= i; idx_row++) {
            solution->searchDir[idx_row - 1] *= rtInf;
          }
        }
        break;
      }
    } else {
      nullStartIdx_tmp = 7 * qrmanager->ncols + 1;
      switch (objective->objtype) {
       case 5:
        for (order = 0; order < mNull_tmp; order++) {
          memspace->workspace_double[order] = -qrmanager->Q[nVar + 7 *
            (qrmanager->ncols + order)];
        }

        order = qrmanager->mrows - 1;
        if (qrmanager->mrows != 0) {
          if (0 <= order) {
            memset(&solution->searchDir[0], 0, (order + 1) * sizeof(double));
          }

          ix = 0;
          i = nullStartIdx_tmp + 7 * (mNull_tmp - 1);
          for (iac = nullStartIdx_tmp; iac <= i; iac += 7) {
            nVars = 0;
            idx_row = iac + order;
            for (ia = iac; ia <= idx_row; ia++) {
              solution->searchDir[nVars] += qrmanager->Q[ia - 1] *
                memspace->workspace_double[ix];
              nVars++;
            }

            ix++;
          }
        }
        break;

       default:
        switch (objective->objtype) {
         case 3:
          xgemm(qrmanager->mrows, mNull_tmp, qrmanager->mrows, H,
                qrmanager->mrows, qrmanager->Q, nullStartIdx_tmp,
                memspace->workspace_double);
          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q,
                  nullStartIdx_tmp, memspace->workspace_double,
                  cholmanager->FMat);
          break;

         default:
          nVars = qrmanager->mrows;
          xgemm(objective->nvar, mNull_tmp, objective->nvar, H, objective->nvar,
                qrmanager->Q, nullStartIdx_tmp, memspace->workspace_double);
          for (order = 0; order < mNull_tmp; order++) {
            i = objective->nvar + 1;
            for (idx_row = i; idx_row <= nVars; idx_row++) {
              memspace->workspace_double[(idx_row + 7 * order) - 1] = 0.0 *
                qrmanager->Q[(idx_row + 7 * (order + qrmanager->ncols)) - 1];
            }
          }

          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q,
                  nullStartIdx_tmp, memspace->workspace_double,
                  cholmanager->FMat);
          break;
        }

        SCALED_REG_PRIMAL = 1.4901161193847656E-6 * (double)mNull_tmp;
        cholmanager->ndims = mNull_tmp;
        nVars = 1;
        if (mNull_tmp > 1) {
          ix = 0;
          smax = fabs(cholmanager->FMat[0]);
          for (idx_row = 2; idx_row <= mNull_tmp; idx_row++) {
            ix += 8;
            s = fabs(cholmanager->FMat[ix]);
            if (s > smax) {
              nVars = idx_row;
              smax = s;
            }
          }
        }

        smax = fabs(cholmanager->FMat[(nVars << 3) - 1]) *
          2.2204460492503131E-16;
        if (!(smax > SCALED_REG_PRIMAL)) {
          smax = SCALED_REG_PRIMAL;
        }

        cholmanager->regTol_ = smax;
        if (mNull_tmp > 128) {
          idx_row = 0;
          exitg1 = false;
          while ((!exitg1) && (idx_row < mNull_tmp)) {
            nVars = (idx_row << 3) + 1;
            order = mNull_tmp - idx_row;
            if (idx_row + 48 <= mNull_tmp) {
              partialColLDL3_(cholmanager, nVars, order, SCALED_REG_PRIMAL);
              idx_row += 48;
            } else {
              fullColLDL2_(cholmanager, nVars, order, SCALED_REG_PRIMAL);
              exitg1 = true;
            }
          }
        } else {
          fullColLDL2_(cholmanager, 1, mNull_tmp, SCALED_REG_PRIMAL);
        }

        if (cholmanager->ConvexCheck) {
          order = 0;
          do {
            exitg2 = 0;
            if (order <= mNull_tmp - 1) {
              if (cholmanager->FMat[order + 7 * order] <= 0.0) {
                cholmanager->info = -order - 1;
                exitg2 = 1;
              } else {
                order++;
              }
            } else {
              cholmanager->ConvexCheck = false;
              exitg2 = 1;
            }
          } while (exitg2 == 0);
        }

        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          if (qrmanager->mrows != 0) {
            memset(&memspace->workspace_double[0], 0, mNull_tmp * sizeof(double));
            nVars = 0;
            i = nullStartIdx_tmp + 7 * (mNull_tmp - 1);
            for (iac = nullStartIdx_tmp; iac <= i; iac += 7) {
              ix = 0;
              smax = 0.0;
              idx_row = iac + nVar;
              for (ia = iac; ia <= idx_row; ia++) {
                smax += qrmanager->Q[ia - 1] * objective->grad[ix];
                ix++;
              }

              memspace->workspace_double[nVars] += -smax;
              nVars++;
            }
          }

          nVars = cholmanager->ndims - 2;
          if (cholmanager->ndims != 0) {
            for (idx_row = 0; idx_row <= nVars + 1; idx_row++) {
              order = idx_row + idx_row * 7;
              i = nVars - idx_row;
              for (iac = 0; iac <= i; iac++) {
                ix = (idx_row + iac) + 1;
                memspace->workspace_double[ix] -= memspace->
                  workspace_double[idx_row] * cholmanager->FMat[(order + iac) +
                  1];
              }
            }
          }

          i = cholmanager->ndims;
          for (order = 0; order < i; order++) {
            memspace->workspace_double[order] /= cholmanager->FMat[order + 7 *
              order];
          }

          nVars = cholmanager->ndims;
          if (cholmanager->ndims != 0) {
            for (idx_row = nVars; idx_row >= 1; idx_row--) {
              order = (idx_row - 1) * 7;
              smax = memspace->workspace_double[idx_row - 1];
              i = idx_row + 1;
              for (iac = nVars; iac >= i; iac--) {
                smax -= cholmanager->FMat[(order + iac) - 1] *
                  memspace->workspace_double[iac - 1];
              }

              memspace->workspace_double[idx_row - 1] = smax;
            }
          }

          order = qrmanager->mrows - 1;
          if (qrmanager->mrows != 0) {
            if (0 <= order) {
              memset(&solution->searchDir[0], 0, (order + 1) * sizeof(double));
            }

            ix = 0;
            i = nullStartIdx_tmp + 7 * (mNull_tmp - 1);
            for (iac = nullStartIdx_tmp; iac <= i; iac += 7) {
              nVars = 0;
              idx_row = iac + order;
              for (ia = iac; ia <= idx_row; ia++) {
                solution->searchDir[nVars] += qrmanager->Q[ia - 1] *
                  memspace->workspace_double[ix];
                nVars++;
              }

              ix++;
            }
          }
        }
        break;
      }
    }
  }
}

/*
 * Arguments    : double workspace[28]
 *                struct_T *solution
 *                const g_struct_T *objective
 *                const e_struct_T *qrmanager
 * Return Type  : void
 */
static void compute_lambda(double workspace[28], struct_T *solution, const
  g_struct_T *objective, const e_struct_T *qrmanager)
{
  double c;
  int i;
  int ia;
  int iac;
  int idx;
  int ix;
  int j;
  int nActiveConstr;
  boolean_T guard1 = false;
  boolean_T nonDegenerate;
  nActiveConstr = qrmanager->ncols;
  if (qrmanager->ncols > 0) {
    c = 100.0 * (double)qrmanager->mrows * 2.2204460492503131E-16;
    if ((qrmanager->mrows > 0) && (qrmanager->ncols > 0)) {
      nonDegenerate = true;
    } else {
      nonDegenerate = false;
    }

    if (nonDegenerate) {
      idx = qrmanager->ncols;
      guard1 = false;
      if (qrmanager->mrows < qrmanager->ncols) {
        while ((idx > qrmanager->mrows) && (fabs(qrmanager->QR[(qrmanager->mrows
                  + 7 * (idx - 1)) - 1]) >= c)) {
          idx--;
        }

        nonDegenerate = (idx == qrmanager->mrows);
        if (nonDegenerate) {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        while ((idx >= 1) && (fabs(qrmanager->QR[(idx + 7 * (idx - 1)) - 1]) >=
                              c)) {
          idx--;
        }

        nonDegenerate = (idx == 0);
      }
    }

    if (!nonDegenerate) {
      solution->state = -7;
    } else {
      if (qrmanager->mrows != 0) {
        idx = qrmanager->ncols;
        if (0 <= idx - 1) {
          memset(&workspace[0], 0, idx * sizeof(double));
        }

        idx = 0;
        j = 7 * (qrmanager->ncols - 1) + 1;
        for (iac = 1; iac <= j; iac += 7) {
          ix = 0;
          c = 0.0;
          i = (iac + qrmanager->mrows) - 1;
          for (ia = iac; ia <= i; ia++) {
            c += qrmanager->Q[ia - 1] * objective->grad[ix];
            ix++;
          }

          workspace[idx] += c;
          idx++;
        }
      }

      for (j = nActiveConstr; j >= 1; j--) {
        idx = (j + (j - 1) * 7) - 1;
        workspace[j - 1] /= qrmanager->QR[idx];
        for (i = 0; i <= j - 2; i++) {
          ix = (j - i) - 2;
          workspace[ix] -= workspace[j - 1] * qrmanager->QR[(idx - i) - 1];
        }
      }

      for (idx = 0; idx < nActiveConstr; idx++) {
        solution->lambda[idx] = -workspace[idx];
      }
    }
  }
}

/*
 * Arguments    : int x[7]
 *                int xLen
 *                int workspace[7]
 *                int xMin
 *                int xMax
 * Return Type  : void
 */
static void countsort(int x[7], int xLen, int workspace[7], int xMin, int xMax)
{
  int idx;
  int idxEnd;
  int idxFill;
  int idxStart;
  int maxOffset;
  if ((xLen > 1) && (xMax > xMin)) {
    idxStart = xMax - xMin;
    if (0 <= idxStart) {
      memset(&workspace[0], 0, (idxStart + 1) * sizeof(int));
    }

    maxOffset = idxStart - 1;
    for (idx = 0; idx < xLen; idx++) {
      idxStart = x[idx] - xMin;
      workspace[idxStart]++;
    }

    for (idx = 2; idx <= maxOffset + 2; idx++) {
      workspace[idx - 1] += workspace[idx - 2];
    }

    idxStart = 1;
    idxEnd = workspace[0];
    for (idx = 0; idx <= maxOffset; idx++) {
      for (idxFill = idxStart; idxFill <= idxEnd; idxFill++) {
        x[idxFill - 1] = idx + xMin;
      }

      idxStart = workspace[idx] + 1;
      idxEnd = workspace[idx + 1];
    }

    for (idx = idxStart; idx <= idxEnd; idx++) {
      x[idx - 1] = xMax;
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                int idx
 * Return Type  : void
 */
static void deleteColMoveEnd(e_struct_T *obj, int idx)
{
  double x[49];
  double c;
  double s;
  double temp;
  int b_i;
  int b_ix;
  int b_k;
  int endIdx;
  int i;
  int ix;
  int k;
  int n;
  if (obj->usedPivoting) {
    i = 1;
    while ((i <= obj->ncols) && (obj->jpvt[i - 1] != idx)) {
      i++;
    }

    idx = i;
  }

  if (idx >= obj->ncols) {
    obj->ncols--;
  } else {
    obj->jpvt[idx - 1] = obj->jpvt[obj->ncols - 1];
    b_i = obj->minRowCol;
    for (k = 0; k < b_i; k++) {
      obj->QR[k + 7 * (idx - 1)] = obj->QR[k + 7 * (obj->ncols - 1)];
    }

    obj->ncols--;
    ix = obj->mrows;
    i = obj->ncols;
    if (ix < i) {
      i = ix;
    }

    obj->minRowCol = i;
    if (idx < obj->mrows) {
      ix = obj->mrows - 1;
      endIdx = obj->ncols;
      if (ix < endIdx) {
        endIdx = ix;
      }

      for (k = endIdx; k >= idx; k--) {
        b_i = k + 7 * (idx - 1);
        temp = obj->QR[b_i];
        xrotg(&obj->QR[(k + 7 * (idx - 1)) - 1], &temp, &c, &s);
        obj->QR[b_i] = temp;
        b_ix = 7 * (k - 1);
        obj->QR[k + b_ix] = 0.0;
        i = k + 7 * idx;
        n = obj->ncols - idx;
        memcpy(&x[0], &obj->QR[0], 49U * sizeof(double));
        if (n >= 1) {
          ix = i - 1;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * x[ix] + s * x[i];
            x[i] = c * x[i] - s * x[ix];
            x[ix] = temp;
            i += 7;
            ix += 7;
          }
        }

        n = obj->mrows;
        for (b_i = 0; b_i < 49; b_i++) {
          obj->QR[b_i] = x[b_i];
          x[b_i] = obj->Q[b_i];
        }

        if (obj->mrows >= 1) {
          i = b_ix + 7;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * x[b_ix] + s * x[i];
            x[i] = c * x[i] - s * x[b_ix];
            x[b_ix] = temp;
            i++;
            b_ix++;
          }
        }

        memcpy(&obj->Q[0], &x[0], 49U * sizeof(double));
      }

      b_i = idx + 1;
      for (k = b_i; k <= endIdx; k++) {
        b_ix = 7 * (k - 1);
        i = k + b_ix;
        temp = obj->QR[i];
        xrotg(&obj->QR[(k + 7 * (k - 1)) - 1], &temp, &c, &s);
        obj->QR[i] = temp;
        i = k << 3;
        n = obj->ncols - k;
        memcpy(&x[0], &obj->QR[0], 49U * sizeof(double));
        if (n >= 1) {
          ix = i - 1;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * x[ix] + s * x[i];
            x[i] = c * x[i] - s * x[ix];
            x[ix] = temp;
            i += 7;
            ix += 7;
          }
        }

        n = obj->mrows;
        for (i = 0; i < 49; i++) {
          obj->QR[i] = x[i];
          x[i] = obj->Q[i];
        }

        if (obj->mrows >= 1) {
          i = b_ix + 7;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * x[b_ix] + s * x[i];
            x[i] = c * x[i] - s * x[b_ix];
            x[b_ix] = temp;
            i++;
            b_ix++;
          }
        }

        memcpy(&obj->Q[0], &x[0], 49U * sizeof(double));
      }
    }
  }
}

/*
 * Arguments    : const double H[9]
 *                const double f[3]
 *                struct_T *solution
 *                c_struct_T *memspace
 *                d_struct_T *workingset
 *                b_struct_T runTimeOptions
 *                e_struct_T *qrmanager
 *                f_struct_T *cholmanager
 *                g_struct_T *objective
 * Return Type  : void
 */
static void driver(const double H[9], const double f[3], struct_T *solution,
                   c_struct_T *memspace, d_struct_T *workingset, b_struct_T
                   runTimeOptions, e_struct_T *qrmanager, f_struct_T
                   *cholmanager, g_struct_T *objective)
{
  static const char cv[128] = { '\x00', '\x01', '\x02', '\x03', '\x04', '\x05',
    '\x06', '\x07', '\x08', '	', '\x0a', '\x0b', '\x0c', '\x0d', '\x0e', '\x0f',
    '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18',
    '\x19', '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', ' ', '!', '\"', '#',
    '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2',
    '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'a',
    'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
    'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '[', '\\', ']', '^', '_',
    '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
    'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}',
    '~', '\x7f' };

  static const char t3_Algorithm[10] = { 'a', 'c', 't', 'i', 'v', 'e', '-', 's',
    'e', 't' };

  static const char cv1[8] = { 'q', 'u', 'a', 'd', 'p', 'r', 'o', 'g' };

  static const char t3_SolverName[8] = { 'q', 'u', 'a', 'd', 'p', 'r', 'o', 'g'
  };

  static const char t3_FiniteDifferenceType[7] = { 'f', 'o', 'r', 'w', 'a', 'r',
    'd' };

  static const char t3_Display[5] = { 'f', 'i', 'n', 'a', 'l' };

  h_struct_T options;
  double qtb;
  double tol;
  int PHASEONE;
  int b_nVar;
  int exitg2;
  int i;
  int idx;
  int mWorkingFixed;
  int nDepInd;
  int nVar;
  boolean_T exitg1;
  boolean_T guard1 = false;
  boolean_T okWorkingSet;
  objective->grad[0] = 0.0;
  objective->grad[1] = 0.0;
  objective->grad[2] = 0.0;
  objective->grad[3] = 0.0;
  objective->Hx[0] = 0.0;
  objective->Hx[1] = 0.0;
  objective->Hx[2] = 0.0;
  objective->hasLinear = true;
  objective->nvar = 3;
  objective->maxVar = 4;
  objective->beta = 0.0;
  objective->rho = 0.0;
  objective->objtype = 3;
  objective->prev_objtype = 3;
  objective->prev_nvar = 0;
  objective->prev_hasLinear = false;
  objective->gammaScalar = 0.0;
  memset(&cholmanager->FMat[0], 0, 49U * sizeof(double));
  cholmanager->ldm = 7;
  cholmanager->ndims = 0;
  cholmanager->info = 0;
  cholmanager->scaleFactor = 100.0;
  cholmanager->ConvexCheck = true;
  cholmanager->regTol_ = 0.0;
  memset(&cholmanager->workspace_[0], 0, 336U * sizeof(double));
  memset(&cholmanager->workspace2_[0], 0, 336U * sizeof(double));
  solution->iterations = 0;
  runTimeOptions.RemainFeasible = true;
  nVar = workingset->nVar - 1;
  i = workingset->sizes[0];
  for (idx = 0; idx < i; idx++) {
    solution->xstar[workingset->indexFixed[idx] - 1] = workingset->ub
      [workingset->indexFixed[idx] - 1];
  }

  i = workingset->sizes[3];
  for (idx = 0; idx < i; idx++) {
    if (workingset->isActiveConstr[(workingset->isActiveIdx[3] + idx) - 1]) {
      solution->xstar[workingset->indexLB[idx] - 1] = -workingset->lb
        [workingset->indexLB[idx] - 1];
    }
  }

  i = workingset->sizes[4];
  for (idx = 0; idx < i; idx++) {
    if (workingset->isActiveConstr[(workingset->isActiveIdx[4] + idx) - 1]) {
      solution->xstar[workingset->indexUB[idx] - 1] = workingset->ub
        [workingset->indexUB[idx] - 1];
    }
  }

  solution->state = 82;
  qrmanager->ldq = 7;
  memset(&qrmanager->QR[0], 0, 49U * sizeof(double));
  memset(&qrmanager->Q[0], 0, 49U * sizeof(double));
  qrmanager->mrows = 0;
  qrmanager->ncols = 0;
  for (i = 0; i < 7; i++) {
    qrmanager->jpvt[i] = 0;
    qrmanager->tau[i] = 0.0;
  }

  qrmanager->minRowCol = 0;
  qrmanager->usedPivoting = false;
  b_nVar = workingset->nVar - 1;
  mWorkingFixed = workingset->nWConstr[0] - 1;
  nDepInd = 0;
  if (workingset->nWConstr[0] > 0) {
    for (i = 0; i <= mWorkingFixed; i++) {
      for (PHASEONE = 0; PHASEONE <= b_nVar; PHASEONE++) {
        qrmanager->QR[i + 7 * PHASEONE] = workingset->ATwset[PHASEONE + (i << 2)];
      }
    }

    nDepInd = 0;
    if (0 <= b_nVar) {
      memset(&qrmanager->jpvt[0], 0, (b_nVar + 1) * sizeof(int));
    }

    qrmanager->usedPivoting = true;
    qrmanager->mrows = workingset->nWConstr[0];
    qrmanager->ncols = workingset->nVar;
    b_nVar = workingset->nWConstr[0];
    i = workingset->nVar;
    if (b_nVar < i) {
      i = b_nVar;
    }

    qrmanager->minRowCol = i;
    xgeqp3(qrmanager->QR, workingset->nWConstr[0], workingset->nVar,
           qrmanager->jpvt, qrmanager->tau);
    tol = 100.0 * (double)workingset->nVar * 2.2204460492503131E-16;
    b_nVar = workingset->nVar;
    i = workingset->nWConstr[0];
    if (b_nVar < i) {
      i = b_nVar;
    }

    while ((i > 0) && (fabs(qrmanager->QR[(i + 7 * (i - 1)) - 1]) < tol)) {
      i--;
      nDepInd++;
    }

    if (nDepInd > 0) {
      computeQ_(qrmanager, workingset->nWConstr[0]);
      idx = 0;
      exitg1 = false;
      while ((!exitg1) && (idx <= nDepInd - 1)) {
        i = 7 * (mWorkingFixed - idx);
        qtb = 0.0;
        b_nVar = 0;
        for (PHASEONE = 0; PHASEONE <= mWorkingFixed; PHASEONE++) {
          qtb += qrmanager->Q[i] * workingset->bwset[b_nVar];
          i++;
          b_nVar++;
        }

        if (fabs(qtb) >= tol) {
          nDepInd = -1;
          exitg1 = true;
        } else {
          idx++;
        }
      }
    }

    if (nDepInd > 0) {
      for (idx = 0; idx <= mWorkingFixed; idx++) {
        qrmanager->jpvt[idx] = 1;
      }

      i = workingset->nWConstr[0] + 1;
      if (i <= mWorkingFixed + 1) {
        memset(&qrmanager->jpvt[i + -1], 0, ((mWorkingFixed - i) + 2) * sizeof
               (int));
      }

      factorQRE(qrmanager, workingset->ATwset, workingset->nVar,
                workingset->nWConstr[0]);
      for (idx = 0; idx < nDepInd; idx++) {
        memspace->workspace_int[idx] = qrmanager->jpvt[((mWorkingFixed - nDepInd)
          + idx) + 1];
      }

      countsort(memspace->workspace_int, nDepInd, memspace->workspace_sort, 1,
                workingset->nWConstr[0]);
      for (idx = nDepInd; idx >= 1; idx--) {
        removeEqConstr(workingset, memspace->workspace_int[idx - 1]);
      }
    }
  }

  if (nDepInd != -1) {
    RemoveDependentIneq_(workingset, qrmanager, memspace, 100.0);
    okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_double,
      solution->xstar, workingset, qrmanager);
    guard1 = false;
    if (!okWorkingSet) {
      RemoveDependentIneq_(workingset, qrmanager, memspace, 1000.0);
      okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_double,
        solution->xstar, workingset, qrmanager);
      if (!okWorkingSet) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 && (workingset->nWConstr[0] + workingset->nWConstr[1] ==
                   workingset->nVar)) {
      tol = b_maxConstraintViolation(workingset, solution->xstar);
      if (tol > 1.0E-8) {
        solution->state = -2;
      }
    }
  } else {
    solution->state = -3;
    i = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
    b_nVar = workingset->nActiveConstr;
    for (nDepInd = i; nDepInd <= b_nVar; nDepInd++) {
      workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
        Wid[nDepInd - 1] - 1] + workingset->Wlocalidx[nDepInd - 1]) - 2] = false;
    }

    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  }

  options.InitDamping = 0.01;
  for (i = 0; i < 7; i++) {
    options.FiniteDifferenceType[i] = t3_FiniteDifferenceType[i];
  }

  options.SpecifyObjectiveGradient = false;
  options.ScaleProblem = false;
  options.SpecifyConstraintGradient = false;
  options.NonFiniteSupport = true;
  options.IterDisplaySQP = false;
  options.FiniteDifferenceStepSize = -1.0;
  options.MaxFunctionEvaluations = -1.0;
  options.IterDisplayQP = false;
  options.PricingTolerance = 0.0;
  for (i = 0; i < 10; i++) {
    options.Algorithm[i] = t3_Algorithm[i];
  }

  options.ObjectiveLimit = -1.0E+20;
  options.ConstraintTolerance = 1.0E-8;
  options.OptimalityTolerance = 1.0E-8;
  options.StepTolerance = 1.0E-8;
  options.MaxIterations = -1.0;
  options.FunctionTolerance = rtInf;
  for (i = 0; i < 8; i++) {
    options.SolverName[i] = t3_SolverName[i];
  }

  options.CheckGradients = false;
  options.Diagnostics[0] = 'o';
  options.Diagnostics[1] = 'f';
  options.Diagnostics[2] = 'f';
  options.DiffMaxChange = rtInf;
  options.DiffMinChange = 0.0;
  for (i = 0; i < 5; i++) {
    options.Display[i] = t3_Display[i];
  }

  options.FunValCheck[0] = 'o';
  options.FunValCheck[1] = 'f';
  options.FunValCheck[2] = 'f';
  options.UseParallel = false;
  options.LinearSolver[0] = 'a';
  options.LinearSolver[1] = 'u';
  options.LinearSolver[2] = 't';
  options.LinearSolver[3] = 'o';
  options.SubproblemAlgorithm[0] = 'c';
  options.SubproblemAlgorithm[1] = 'g';
  if (solution->state >= 0) {
    solution->iterations = 0;
    solution->maxConstr = b_maxConstraintViolation(workingset, solution->xstar);
    guard1 = false;
    if (solution->maxConstr > 1.0E-8) {
      phaseone(H, f, solution, memspace, workingset, qrmanager, &options,
               &runTimeOptions, cholmanager, objective);
      if (solution->state != 0) {
        solution->maxConstr = b_maxConstraintViolation(workingset,
          solution->xstar);
        if (solution->maxConstr > 1.0E-8) {
          for (PHASEONE = 0; PHASEONE < 7; PHASEONE++) {
            solution->lambda[PHASEONE] = 0.0;
          }

          solution->fstar = computeFval(objective, memspace->workspace_double, H,
            f, solution->xstar);
          solution->state = -2;
        } else {
          if (solution->maxConstr > 0.0) {
            if (0 <= nVar) {
              memcpy(&solution->searchDir[0], &solution->xstar[0], (nVar + 1) *
                     sizeof(double));
            }

            PresolveWorkingSet(solution, memspace, workingset, qrmanager);
            tol = b_maxConstraintViolation(workingset, solution->xstar);
            if (tol >= solution->maxConstr) {
              solution->maxConstr = tol;
              if (0 <= nVar) {
                memcpy(&solution->xstar[0], &solution->searchDir[0], (nVar + 1) *
                       sizeof(double));
              }
            }
          }

          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
              objective, options.ObjectiveLimit, options.StepTolerance,
              runTimeOptions.MaxIterations, runTimeOptions.ProbRelTolFactor,
              true);
      okWorkingSet = false;
      i = 0;
      do {
        exitg2 = 0;
        if (i < 8) {
          if (cv[(unsigned char)options.SolverName[i]] != cv[(int)cv1[i]]) {
            exitg2 = 1;
          } else {
            i++;
          }
        } else {
          okWorkingSet = true;
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (okWorkingSet && (solution->state != -6)) {
        solution->maxConstr = b_maxConstraintViolation(workingset,
          solution->xstar);
        computeFirstOrderOpt(solution, objective, workingset->nVar,
                             workingset->ATwset, workingset->nActiveConstr,
                             memspace->workspace_double);
        while ((solution->iterations < runTimeOptions.MaxIterations) &&
               ((solution->state == -7) || ((solution->state == 1) &&
                 ((solution->maxConstr > 1.0E-8) || (solution->firstorderopt >
                   1.0E-8 * runTimeOptions.ProbRelTolFactor))))) {
          feasibleX0ForWorkingSet(memspace->workspace_double, solution->xstar,
            workingset, qrmanager);
          PresolveWorkingSet(solution, memspace, workingset, qrmanager);
          mWorkingFixed = workingset->probType;
          nVar = workingset->nVar;
          solution->xstar[3] = solution->maxConstr + 1.0;
          if (workingset->probType == 3) {
            PHASEONE = 1;
          } else {
            PHASEONE = 4;
          }

          i = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          b_nVar = workingset->nActiveConstr;
          for (nDepInd = i; nDepInd <= b_nVar; nDepInd++) {
            workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
              Wid[nDepInd - 1] - 1] + workingset->Wlocalidx[nDepInd - 1]) - 2] =
              false;
          }

          workingset->nWConstr[2] = 0;
          workingset->nWConstr[3] = 0;
          workingset->nWConstr[4] = 0;
          workingset->nActiveConstr = workingset->nWConstr[0] +
            workingset->nWConstr[1];
          setProblemType(workingset, PHASEONE);
          objective->prev_objtype = objective->objtype;
          objective->prev_nvar = objective->nvar;
          objective->prev_hasLinear = objective->hasLinear;
          objective->objtype = 5;
          objective->nvar = 4;
          objective->gammaScalar = 1.0;
          objective->hasLinear = true;
          solution->fstar = computeFval(objective, memspace->workspace_double, H,
            f, solution->xstar);
          solution->state = 5;
          iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
                  objective, 1.0E-8, 1.4901161193847657E-10,
                  runTimeOptions.MaxIterations, runTimeOptions.ProbRelTolFactor,
                  false);
          if (workingset->isActiveConstr[(workingset->isActiveIdx[3] +
               workingset->sizes[3]) - 2]) {
            idx = workingset->sizes[0];
            exitg1 = false;
            while ((!exitg1) && (idx + 1 <= workingset->nActiveConstr)) {
              if ((workingset->Wid[idx] == 4) && (workingset->Wlocalidx[idx] ==
                   workingset->sizes[3])) {
                removeConstr(workingset, idx + 1);
                exitg1 = true;
              } else {
                idx++;
              }
            }
          }

          i = workingset->nActiveConstr;
          b_nVar = workingset->sizes[0];
          while ((i > b_nVar) && (i > nVar)) {
            removeConstr(workingset, i);
            i--;
          }

          solution->maxConstr = solution->xstar[3];
          setProblemType(workingset, mWorkingFixed);
          objective->objtype = objective->prev_objtype;
          objective->nvar = objective->prev_nvar;
          objective->hasLinear = objective->prev_hasLinear;
          iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
                  objective, -1.0E+20, 1.0E-8, runTimeOptions.MaxIterations,
                  runTimeOptions.ProbRelTolFactor, false);
          solution->maxConstr = b_maxConstraintViolation(workingset,
            solution->xstar);
          computeFirstOrderOpt(solution, objective, workingset->nVar,
                               workingset->ATwset, workingset->nActiveConstr,
                               memspace->workspace_double);
        }
      }
    }
  }
}

/*
 * Arguments    : f_struct_T *obj
 *                const double A[9]
 *                int ndims
 *                int ldA
 * Return Type  : void
 */
static void factor(f_struct_T *obj, const double A[9], int ndims, int ldA)
{
  double SCALED_REG_PRIMAL;
  double s;
  double smax;
  int A_maxDiag_idx;
  int exitg2;
  int idx;
  int ix;
  int k;
  boolean_T exitg1;
  SCALED_REG_PRIMAL = 1.4901161193847656E-6 * (double)ndims;
  obj->ndims = ndims;
  for (idx = 0; idx < ndims; idx++) {
    A_maxDiag_idx = ldA * idx;
    ix = 7 * idx;
    for (k = 0; k < ndims; k++) {
      obj->FMat[ix + k] = A[A_maxDiag_idx + k];
    }
  }

  if (ndims < 1) {
    A_maxDiag_idx = 0;
  } else {
    A_maxDiag_idx = 1;
    if (ndims > 1) {
      ix = 0;
      smax = fabs(obj->FMat[0]);
      for (k = 2; k <= ndims; k++) {
        ix += 8;
        s = fabs(obj->FMat[ix]);
        if (s > smax) {
          A_maxDiag_idx = k;
          smax = s;
        }
      }
    }
  }

  smax = fabs(obj->FMat[(A_maxDiag_idx << 3) - 1]) * 2.2204460492503131E-16;
  s = fabs(SCALED_REG_PRIMAL);
  if (smax > s) {
    s = smax;
  }

  obj->regTol_ = s;
  if (ndims > 128) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < ndims)) {
      A_maxDiag_idx = (k << 3) + 1;
      ix = ndims - k;
      if (k + 48 <= ndims) {
        partialColLDL3_(obj, A_maxDiag_idx, ix, SCALED_REG_PRIMAL);
        k += 48;
      } else {
        fullColLDL2_(obj, A_maxDiag_idx, ix, SCALED_REG_PRIMAL);
        exitg1 = true;
      }
    }
  } else {
    fullColLDL2_(obj, 1, ndims, SCALED_REG_PRIMAL);
  }

  if (obj->ConvexCheck) {
    idx = 0;
    do {
      exitg2 = 0;
      if (idx <= ndims - 1) {
        if (obj->FMat[idx + 7 * idx] <= 0.0) {
          obj->info = -idx - 1;
          exitg2 = 1;
        } else {
          idx++;
        }
      } else {
        obj->ConvexCheck = false;
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double A[28]
 *                int mrows
 *                int ncols
 * Return Type  : void
 */
static void factorQR(e_struct_T *obj, const double A[28], int mrows, int ncols)
{
  int i;
  int idx;
  int ix0;
  int k;
  boolean_T guard1 = false;
  ix0 = mrows * ncols;
  guard1 = false;
  if (ix0 > 0) {
    for (idx = 0; idx < ncols; idx++) {
      ix0 = idx << 2;
      i = 7 * idx;
      for (k = 0; k < mrows; k++) {
        obj->QR[i + k] = A[ix0 + k];
      }
    }

    guard1 = true;
  } else if (ix0 == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = false;
    obj->mrows = mrows;
    obj->ncols = ncols;
    for (idx = 0; idx < ncols; idx++) {
      obj->jpvt[idx] = idx + 1;
    }

    if (mrows < ncols) {
      ix0 = mrows;
    } else {
      ix0 = ncols;
    }

    obj->minRowCol = ix0;
    for (i = 0; i < 7; i++) {
      obj->tau[i] = 0.0;
    }

    if (ix0 >= 1) {
      for (i = 0; i < 7; i++) {
        obj->tau[i] = 0.0;
      }

      qrf(obj->QR, mrows, ncols, ix0, obj->tau);
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double A[28]
 *                int mrows
 *                int ncols
 * Return Type  : void
 */
static void factorQRE(e_struct_T *obj, const double A[28], int mrows, int ncols)
{
  int idx;
  int ix0;
  int iy0;
  int k;
  boolean_T guard1 = false;
  ix0 = mrows * ncols;
  guard1 = false;
  if (ix0 > 0) {
    for (idx = 0; idx < ncols; idx++) {
      ix0 = idx << 2;
      iy0 = 7 * idx;
      for (k = 0; k < mrows; k++) {
        obj->QR[iy0 + k] = A[ix0 + k];
      }
    }

    guard1 = true;
  } else if (ix0 == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = true;
    obj->mrows = mrows;
    obj->ncols = ncols;
    if (mrows < ncols) {
      ix0 = mrows;
    } else {
      ix0 = ncols;
    }

    obj->minRowCol = ix0;
    xgeqp3(obj->QR, mrows, ncols, obj->jpvt, obj->tau);
  }
}

/*
 * Arguments    : double workspace[28]
 *                double xCurrent[4]
 *                const d_struct_T *workingset
 *                e_struct_T *qrmanager
 * Return Type  : boolean_T
 */
static boolean_T feasibleX0ForWorkingSet(double workspace[28], double xCurrent[4],
  const d_struct_T *workingset, e_struct_T *qrmanager)
{
  double B[28];
  double c;
  double constrViolation_basicX;
  int b_i;
  int exitg1;
  int i;
  int ia;
  int iac;
  int ix;
  int iy;
  int mWConstr;
  int nVar;
  boolean_T nonDegenerateWset;
  mWConstr = workingset->nActiveConstr - 1;
  nVar = workingset->nVar;
  nonDegenerateWset = true;
  if (workingset->nActiveConstr != 0) {
    for (iy = 0; iy <= mWConstr; iy++) {
      workspace[iy] = workingset->bwset[iy];
      workspace[iy + 7] = workingset->bwset[iy];
    }

    if (workingset->nActiveConstr != 0) {
      iy = 0;
      i = ((workingset->nActiveConstr - 1) << 2) + 1;
      for (iac = 1; iac <= i; iac += 4) {
        ix = 0;
        c = 0.0;
        b_i = (iac + nVar) - 1;
        for (ia = iac; ia <= b_i; ia++) {
          c += workingset->ATwset[ia - 1] * xCurrent[ix];
          ix++;
        }

        workspace[iy] += -c;
        iy++;
      }
    }

    if (workingset->nActiveConstr >= workingset->nVar) {
      for (iy = 0; iy < nVar; iy++) {
        for (ix = 0; ix <= mWConstr; ix++) {
          qrmanager->QR[ix + 7 * iy] = workingset->ATwset[iy + (ix << 2)];
        }
      }

      qrmanager->usedPivoting = false;
      qrmanager->mrows = workingset->nActiveConstr;
      qrmanager->ncols = workingset->nVar;
      i = workingset->nVar;
      for (iy = 0; iy < i; iy++) {
        qrmanager->jpvt[iy] = iy + 1;
      }

      iy = workingset->nActiveConstr;
      ix = workingset->nVar;
      if (iy < ix) {
        ix = iy;
      }

      qrmanager->minRowCol = ix;
      for (b_i = 0; b_i < 7; b_i++) {
        qrmanager->tau[b_i] = 0.0;
      }

      if (ix >= 1) {
        for (b_i = 0; b_i < 7; b_i++) {
          qrmanager->tau[b_i] = 0.0;
        }

        qrf(qrmanager->QR, workingset->nActiveConstr, workingset->nVar, ix,
            qrmanager->tau);
      }

      computeQ_(qrmanager, workingset->nActiveConstr);
      memcpy(&B[0], &workspace[0], 28U * sizeof(double));
      if (1 <= nVar) {
        memset(&workspace[0], 0, nVar * sizeof(double));
      }

      i = nVar + 7;
      if (8 <= i) {
        memset(&workspace[7], 0, (i + -7) * sizeof(double));
      }

      iy = -1;
      for (iac = 1; iac <= nVar; iac++) {
        c = 0.0;
        for (ix = 0; ix <= mWConstr; ix++) {
          c += qrmanager->Q[(ix + iy) + 1] * B[ix];
        }

        workspace[iac - 1] += c;
        iy += 7;
      }

      iy = -1;
      i = nVar + 7;
      for (iac = 8; iac <= i; iac++) {
        c = 0.0;
        for (ix = 0; ix <= mWConstr; ix++) {
          c += qrmanager->Q[(ix + iy) + 1] * B[ix + 7];
        }

        workspace[iac - 1] += c;
        iy += 7;
      }

      for (ix = nVar; ix >= 1; ix--) {
        iy = 7 * (ix - 1) - 1;
        c = workspace[ix + -1];
        if (c != 0.0) {
          workspace[ix + -1] = c / qrmanager->QR[ix + iy];
          for (b_i = 0; b_i <= ix - 2; b_i++) {
            workspace[b_i] -= workspace[ix + -1] * qrmanager->QR[(b_i + iy) + 1];
          }
        }
      }

      for (ix = nVar; ix >= 1; ix--) {
        iy = 7 * (ix - 1) - 1;
        c = workspace[ix + 6];
        if (c != 0.0) {
          workspace[ix + 6] = c / qrmanager->QR[ix + iy];
          for (b_i = 0; b_i <= ix - 2; b_i++) {
            workspace[b_i + 7] -= workspace[ix + 6] * qrmanager->QR[(b_i + iy) +
              1];
          }
        }
      }
    } else {
      factorQR(qrmanager, workingset->ATwset, workingset->nVar,
               workingset->nActiveConstr);
      computeQ_(qrmanager, qrmanager->minRowCol);
      for (b_i = 0; b_i <= mWConstr; b_i++) {
        iy = 7 * b_i;
        c = workspace[b_i];
        for (ix = 0; ix < b_i; ix++) {
          c -= qrmanager->QR[ix + iy] * workspace[ix];
        }

        workspace[b_i] = c / qrmanager->QR[b_i + iy];
      }

      for (b_i = 0; b_i <= mWConstr; b_i++) {
        iy = 7 * b_i;
        c = workspace[b_i + 7];
        for (ix = 0; ix < b_i; ix++) {
          c -= qrmanager->QR[ix + iy] * workspace[ix + 7];
        }

        workspace[b_i + 7] = c / qrmanager->QR[b_i + iy];
      }

      memcpy(&B[0], &workspace[0], 28U * sizeof(double));
      if (1 <= nVar) {
        memset(&workspace[0], 0, nVar * sizeof(double));
      }

      i = nVar + 7;
      if (8 <= i) {
        memset(&workspace[7], 0, (i + -7) * sizeof(double));
      }

      iy = -1;
      i = mWConstr + 1;
      for (ix = 1; ix <= i; ix++) {
        ia = iy;
        for (iac = 1; iac <= nVar; iac++) {
          ia++;
          workspace[iac - 1] += B[ix - 1] * qrmanager->Q[ia];
        }

        iy += 7;
      }

      iy = -1;
      i = mWConstr + 8;
      for (ix = 8; ix <= i; ix++) {
        ia = iy;
        b_i = nVar + 7;
        for (iac = 8; iac <= b_i; iac++) {
          ia++;
          workspace[iac - 1] += B[ix - 1] * qrmanager->Q[ia];
        }

        iy += 7;
      }
    }

    iy = 0;
    do {
      exitg1 = 0;
      if (iy <= nVar - 1) {
        if (rtIsInf(workspace[iy]) || rtIsNaN(workspace[iy])) {
          nonDegenerateWset = false;
          exitg1 = 1;
        } else {
          c = workspace[iy + 7];
          if (rtIsInf(c) || rtIsNaN(c)) {
            nonDegenerateWset = false;
            exitg1 = 1;
          } else {
            iy++;
          }
        }
      } else {
        iy = nVar - 1;
        for (ix = 0; ix <= iy; ix++) {
          workspace[ix] += xCurrent[ix];
        }

        c = maxConstraintViolation(workingset, workspace, 1);
        constrViolation_basicX = maxConstraintViolation(workingset, workspace, 8);
        if ((c <= 2.2204460492503131E-16) || (c < constrViolation_basicX)) {
          if (0 <= nVar - 1) {
            memcpy(&xCurrent[0], &workspace[0], nVar * sizeof(double));
          }
        } else {
          if (0 <= nVar - 1) {
            memcpy(&xCurrent[0], &workspace[7], nVar * sizeof(double));
          }
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return nonDegenerateWset;
}

/*
 * Arguments    : const double solution_xstar[4]
 *                const double solution_searchDir[4]
 *                int workingset_nVar
 *                const double workingset_lb[4]
 *                const double workingset_ub[4]
 *                const int workingset_indexLB[4]
 *                const int workingset_indexUB[4]
 *                const int workingset_sizes[5]
 *                const int workingset_isActiveIdx[6]
 *                const boolean_T workingset_isActiveConstr[7]
 *                const int workingset_nWConstr[5]
 *                boolean_T isPhaseOne
 *                double *alpha
 *                boolean_T *newBlocking
 *                int *constrType
 *                int *constrIdx
 * Return Type  : void
 */
static void feasibleratiotest(const double solution_xstar[4], const double
  solution_searchDir[4], int workingset_nVar, const double workingset_lb[4],
  const double workingset_ub[4], const int workingset_indexLB[4], const int
  workingset_indexUB[4], const int workingset_sizes[5], const int
  workingset_isActiveIdx[6], const boolean_T workingset_isActiveConstr[7], const
  int workingset_nWConstr[5], boolean_T isPhaseOne, double *alpha, boolean_T
  *newBlocking, int *constrType, int *constrIdx)
{
  double denomTol;
  double phaseOneCorrectionP;
  double phaseOneCorrectionX;
  double pk_corrected;
  double ratio;
  double u0;
  int i;
  int i1;
  int idx;
  int totalUB;
  totalUB = workingset_sizes[4];
  *alpha = 1.0E+30;
  *newBlocking = false;
  *constrType = 0;
  *constrIdx = 0;
  denomTol = 2.2204460492503131E-13 * b_xnrm2(workingset_nVar,
    solution_searchDir);
  if (workingset_nWConstr[3] < workingset_sizes[3]) {
    phaseOneCorrectionX = (double)isPhaseOne * solution_xstar[workingset_nVar -
      1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir[workingset_nVar - 1];
    i = workingset_sizes[3];
    for (idx = 0; idx <= i - 2; idx++) {
      i1 = workingset_indexLB[idx];
      pk_corrected = -solution_searchDir[i1 - 1] - phaseOneCorrectionP;
      if ((pk_corrected > denomTol) && (!workingset_isActiveConstr
           [(workingset_isActiveIdx[3] + idx) - 1])) {
        ratio = (-solution_xstar[i1 - 1] - workingset_lb[i1 - 1]) -
          phaseOneCorrectionX;
        u0 = fabs(ratio);
        if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
          u0 = 1.0E-8 - ratio;
        }

        pk_corrected = u0 / pk_corrected;
        if (pk_corrected < *alpha) {
          *alpha = pk_corrected;
          *constrType = 4;
          *constrIdx = idx + 1;
          *newBlocking = true;
        }
      }
    }

    i = workingset_indexLB[workingset_sizes[3] - 1] - 1;
    pk_corrected = -solution_searchDir[i];
    if ((pk_corrected > denomTol) && (!workingset_isActiveConstr
         [(workingset_isActiveIdx[3] + workingset_sizes[3]) - 2])) {
      ratio = -solution_xstar[i] - workingset_lb[i];
      u0 = fabs(ratio);
      if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
        u0 = 1.0E-8 - ratio;
      }

      pk_corrected = u0 / pk_corrected;
      if (pk_corrected < *alpha) {
        *alpha = pk_corrected;
        *constrType = 4;
        *constrIdx = workingset_sizes[3];
        *newBlocking = true;
      }
    }
  }

  if (workingset_nWConstr[4] < workingset_sizes[4]) {
    phaseOneCorrectionX = (double)isPhaseOne * solution_xstar[workingset_nVar -
      1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir[workingset_nVar - 1];
    for (idx = 0; idx < totalUB; idx++) {
      i = workingset_indexUB[idx];
      pk_corrected = solution_searchDir[i - 1] - phaseOneCorrectionP;
      if ((pk_corrected > denomTol) && (!workingset_isActiveConstr
           [(workingset_isActiveIdx[4] + idx) - 1])) {
        ratio = (solution_xstar[i - 1] - workingset_ub[i - 1]) -
          phaseOneCorrectionX;
        u0 = fabs(ratio);
        if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
          u0 = 1.0E-8 - ratio;
        }

        pk_corrected = u0 / pk_corrected;
        if (pk_corrected < *alpha) {
          *alpha = pk_corrected;
          *constrType = 5;
          *constrIdx = idx + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (!isPhaseOne) {
    if ((*newBlocking) && (*alpha > 1.0)) {
      *newBlocking = false;
    }

    if (!(*alpha < 1.0)) {
      *alpha = 1.0;
    }
  }
}

/*
 * Arguments    : f_struct_T *obj
 *                int LD_offset
 *                int NColsRemain
 *                double REG_PRIMAL
 * Return Type  : void
 */
static void fullColLDL2_(f_struct_T *obj, int LD_offset, int NColsRemain, double
  REG_PRIMAL)
{
  double alpha1;
  double temp;
  int LD_diagOffset;
  int i;
  int i1;
  int ijA;
  int ix;
  int j;
  int jA;
  int jy;
  int k;
  int offset1;
  int subMatrixDim;
  for (k = 0; k < NColsRemain; k++) {
    LD_diagOffset = (LD_offset + (k << 3)) - 1;
    if (fabs(obj->FMat[LD_diagOffset]) <= obj->regTol_) {
      obj->FMat[LD_diagOffset] += REG_PRIMAL;
    }

    alpha1 = -1.0 / obj->FMat[LD_diagOffset];
    subMatrixDim = (NColsRemain - k) - 2;
    offset1 = LD_diagOffset + 2;
    for (jA = 0; jA <= subMatrixDim; jA++) {
      obj->workspace_[jA] = obj->FMat[(LD_diagOffset + jA) + 1];
    }

    if (!(alpha1 == 0.0)) {
      jA = LD_diagOffset;
      jy = 0;
      for (j = 0; j <= subMatrixDim; j++) {
        if (obj->workspace_[jy] != 0.0) {
          temp = obj->workspace_[jy] * alpha1;
          ix = 0;
          i = jA + 9;
          i1 = subMatrixDim + jA;
          for (ijA = i; ijA <= i1 + 9; ijA++) {
            obj->FMat[ijA - 1] += obj->workspace_[ix] * temp;
            ix++;
          }
        }

        jy++;
        jA += 7;
      }
    }

    alpha1 = 1.0 / obj->FMat[LD_diagOffset];
    i = LD_diagOffset + subMatrixDim;
    for (jA = offset1; jA <= i + 2; jA++) {
      obj->FMat[jA - 1] *= alpha1;
    }
  }

  jA = (LD_offset + ((NColsRemain - 1) << 3)) - 1;
  if (fabs(obj->FMat[jA]) <= obj->regTol_) {
    obj->FMat[jA] += REG_PRIMAL;
  }
}

/*
 * Arguments    : const double H[9]
 *                const double f[3]
 *                struct_T *solution
 *                c_struct_T *memspace
 *                d_struct_T *workingset
 *                e_struct_T *qrmanager
 *                f_struct_T *cholmanager
 *                g_struct_T *objective
 *                double options_ObjectiveLimit
 *                double options_StepTolerance
 *                int runTimeOptions_MaxIterations
 *                double runTimeOptions_ProbRelTolFactor
 *                boolean_T runTimeOptions_RemainFeasible
 * Return Type  : void
 */
static void iterate(const double H[9], const double f[3], struct_T *solution,
                    c_struct_T *memspace, d_struct_T *workingset, e_struct_T
                    *qrmanager, f_struct_T *cholmanager, g_struct_T *objective,
                    double options_ObjectiveLimit, double options_StepTolerance,
                    int runTimeOptions_MaxIterations, double
                    runTimeOptions_ProbRelTolFactor, boolean_T
                    runTimeOptions_RemainFeasible)
{
  double alpha;
  double denomTol;
  double p_max;
  double phaseOneCorrectionP;
  double phaseOneCorrectionX;
  double pk_corrected;
  double ratio;
  double ratio_tmp;
  double tolDelta;
  double u0;
  int TYPE;
  int activeConstrChangedType;
  int activeSetChangeID;
  int exitg1;
  int globalActiveConstrIdx;
  int i;
  int idx;
  int ixlast;
  int localActiveConstrIdx;
  int nVar;
  int totalUB;
  boolean_T guard1 = false;
  boolean_T newBlocking;
  boolean_T subProblemChanged;
  boolean_T updateFval;
  subProblemChanged = true;
  updateFval = true;
  activeSetChangeID = 0;
  TYPE = objective->objtype;
  tolDelta = 6.7434957617430445E-7;
  nVar = workingset->nVar;
  globalActiveConstrIdx = 0;
  computeGrad_StoreHx(objective, H, f, solution->xstar);
  solution->fstar = computeFval_ReuseHx(objective, memspace->workspace_double, f,
    solution->xstar);
  if (solution->iterations < runTimeOptions_MaxIterations) {
    solution->state = -5;
  } else {
    solution->state = 0;
  }

  for (totalUB = 0; totalUB < 7; totalUB++) {
    solution->lambda[totalUB] = 0.0;
  }

  do {
    exitg1 = 0;
    if (solution->state == -5) {
      guard1 = false;
      if (subProblemChanged) {
        switch (activeSetChangeID) {
         case 1:
          squareQ_appendCol(qrmanager, workingset->ATwset,
                            ((workingset->nActiveConstr - 1) << 2) + 1);
          break;

         case -1:
          deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
          break;

         default:
          factorQR(qrmanager, workingset->ATwset, nVar,
                   workingset->nActiveConstr);
          computeQ_(qrmanager, qrmanager->mrows);
          break;
        }

        compute_deltax(H, solution, memspace, qrmanager, cholmanager, objective);
        if (solution->state != -5) {
          exitg1 = 1;
        } else if ((b_xnrm2(nVar, solution->searchDir) < options_StepTolerance) ||
                   (workingset->nActiveConstr >= nVar)) {
          guard1 = true;
        } else {
          updateFval = (TYPE == 5);
          if (updateFval || runTimeOptions_RemainFeasible) {
            feasibleratiotest(solution->xstar, solution->searchDir,
                              workingset->nVar, workingset->lb, workingset->ub,
                              workingset->indexLB, workingset->indexUB,
                              workingset->sizes, workingset->isActiveIdx,
                              workingset->isActiveConstr, workingset->nWConstr,
                              updateFval, &alpha, &newBlocking,
                              &activeConstrChangedType, &localActiveConstrIdx);
          } else {
            totalUB = workingset->sizes[4];
            alpha = 1.0E+30;
            newBlocking = false;
            activeConstrChangedType = 0;
            localActiveConstrIdx = 0;
            p_max = 0.0;
            denomTol = 2.2204460492503131E-13 * b_xnrm2(workingset->nVar,
              solution->searchDir);
            if (workingset->nWConstr[3] < workingset->sizes[3]) {
              ixlast = workingset->nVar - 1;
              phaseOneCorrectionX = 0.0 * solution->xstar[ixlast];
              phaseOneCorrectionP = 0.0 * solution->searchDir[ixlast];
              i = workingset->sizes[3];
              for (idx = 0; idx <= i - 2; idx++) {
                pk_corrected = -solution->searchDir[workingset->indexLB[idx] - 1]
                  - phaseOneCorrectionP;
                if ((pk_corrected > denomTol) && (!workingset->isActiveConstr
                     [(workingset->isActiveIdx[3] + idx) - 1])) {
                  ratio_tmp = -solution->xstar[workingset->indexLB[idx] - 1] -
                    workingset->lb[workingset->indexLB[idx] - 1];
                  ratio = (ratio_tmp - tolDelta) - phaseOneCorrectionX;
                  u0 = fabs(ratio);
                  if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
                    u0 = 1.0E-8 - ratio;
                  }

                  ratio = u0 / pk_corrected;
                  if ((ratio <= alpha) && (fabs(pk_corrected) > p_max)) {
                    alpha = ratio;
                    activeConstrChangedType = 4;
                    localActiveConstrIdx = idx + 1;
                    newBlocking = true;
                  }

                  ratio = ratio_tmp - phaseOneCorrectionX;
                  u0 = fabs(ratio);
                  if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
                    u0 = 1.0E-8 - ratio;
                  }

                  ratio = u0 / pk_corrected;
                  if (ratio < alpha) {
                    alpha = ratio;
                    activeConstrChangedType = 4;
                    localActiveConstrIdx = idx + 1;
                    newBlocking = true;
                    p_max = fabs(pk_corrected);
                  }
                }
              }

              i = workingset->indexLB[workingset->sizes[3] - 1] - 1;
              phaseOneCorrectionX = solution->searchDir[i];
              if ((-phaseOneCorrectionX > denomTol) &&
                  (!workingset->isActiveConstr[(workingset->isActiveIdx[3] +
                    workingset->sizes[3]) - 2])) {
                ratio_tmp = -solution->xstar[i] - workingset->lb[i];
                ratio = ratio_tmp - tolDelta;
                u0 = fabs(ratio);
                if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
                  u0 = 1.0E-8 - ratio;
                }

                ratio = u0 / -phaseOneCorrectionX;
                if ((ratio <= alpha) && (fabs(phaseOneCorrectionX) > p_max)) {
                  alpha = ratio;
                  activeConstrChangedType = 4;
                  localActiveConstrIdx = workingset->sizes[3];
                  newBlocking = true;
                }

                u0 = fabs(ratio_tmp);
                if ((!(u0 < 1.0E-8 - ratio_tmp)) && (!rtIsNaN(1.0E-8 - ratio_tmp)))
                {
                  u0 = 1.0E-8 - ratio_tmp;
                }

                ratio = u0 / -phaseOneCorrectionX;
                if (ratio < alpha) {
                  alpha = ratio;
                  activeConstrChangedType = 4;
                  localActiveConstrIdx = workingset->sizes[3];
                  newBlocking = true;
                  p_max = fabs(solution->searchDir[i]);
                }
              }
            }

            if (workingset->nWConstr[4] < workingset->sizes[4]) {
              ixlast = workingset->nVar - 1;
              phaseOneCorrectionX = 0.0 * solution->xstar[ixlast];
              phaseOneCorrectionP = 0.0 * solution->searchDir[ixlast];
              for (idx = 0; idx < totalUB; idx++) {
                pk_corrected = solution->searchDir[workingset->indexUB[idx] - 1]
                  - phaseOneCorrectionP;
                if ((pk_corrected > denomTol) && (!workingset->isActiveConstr
                     [(workingset->isActiveIdx[4] + idx) - 1])) {
                  ratio_tmp = solution->xstar[workingset->indexUB[idx] - 1] -
                    workingset->ub[workingset->indexUB[idx] - 1];
                  ratio = (ratio_tmp - tolDelta) - phaseOneCorrectionX;
                  u0 = fabs(ratio);
                  if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
                    u0 = 1.0E-8 - ratio;
                  }

                  ratio = u0 / pk_corrected;
                  if ((ratio <= alpha) && (fabs(pk_corrected) > p_max)) {
                    alpha = ratio;
                    activeConstrChangedType = 5;
                    localActiveConstrIdx = idx + 1;
                    newBlocking = true;
                  }

                  ratio = ratio_tmp - phaseOneCorrectionX;
                  u0 = fabs(ratio);
                  if ((!(u0 < 1.0E-8 - ratio)) && (!rtIsNaN(1.0E-8 - ratio))) {
                    u0 = 1.0E-8 - ratio;
                  }

                  ratio = u0 / pk_corrected;
                  if (ratio < alpha) {
                    alpha = ratio;
                    activeConstrChangedType = 5;
                    localActiveConstrIdx = idx + 1;
                    newBlocking = true;
                    p_max = fabs(pk_corrected);
                  }
                }
              }
            }

            tolDelta += 6.608625846508183E-7;
            if (p_max > 0.0) {
              ratio = 6.608625846508183E-7 / p_max;
              if (!(alpha > ratio)) {
                alpha = ratio;
              }
            }

            if (newBlocking && (alpha > 1.0)) {
              newBlocking = false;
            }

            if (!(alpha < 1.0)) {
              alpha = 1.0;
            }
          }

          if (newBlocking) {
            switch (activeConstrChangedType) {
             case 3:
              workingset->nWConstr[2]++;
              workingset->isActiveConstr[(workingset->isActiveIdx[2] +
                localActiveConstrIdx) - 2] = true;
              workingset->nActiveConstr++;
              workingset->Wid[workingset->nActiveConstr - 1] = 3;
              workingset->Wlocalidx[workingset->nActiveConstr - 1] =
                localActiveConstrIdx;

              /* A check that is always false is detected at compile-time. Eliminating code that follows. */
              break;

             case 4:
              addBoundToActiveSetMatrix_(workingset, 4, localActiveConstrIdx);
              break;

             default:
              addBoundToActiveSetMatrix_(workingset, 5, localActiveConstrIdx);
              break;
            }

            activeSetChangeID = 1;
          } else {
            if (objective->objtype == 5) {
              if (b_xnrm2(objective->nvar, solution->searchDir) > 100.0 *
                  (double)objective->nvar * 1.4901161193847656E-8) {
                solution->state = 3;
              } else {
                solution->state = 4;
              }
            }

            subProblemChanged = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }

          if (!(alpha == 0.0)) {
            ixlast = nVar - 1;
            for (totalUB = 0; totalUB <= ixlast; totalUB++) {
              solution->xstar[totalUB] += alpha * solution->searchDir[totalUB];
            }
          }

          computeGrad_StoreHx(objective, H, f, solution->xstar);
          updateFval = true;
          checkStoppingAndUpdateFval(&activeSetChangeID, f, solution, memspace,
            objective, workingset, qrmanager, options_ObjectiveLimit,
            runTimeOptions_MaxIterations, updateFval);
        }
      } else {
        if (0 <= nVar - 1) {
          memset(&solution->searchDir[0], 0, nVar * sizeof(double));
        }

        guard1 = true;
      }

      if (guard1) {
        compute_lambda(memspace->workspace_double, solution, objective,
                       qrmanager);
        if (solution->state != -7) {
          ixlast = 0;
          ratio = 0.0 * runTimeOptions_ProbRelTolFactor * (double)(TYPE != 5);
          i = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          totalUB = workingset->nActiveConstr;
          for (idx = i; idx <= totalUB; idx++) {
            phaseOneCorrectionX = solution->lambda[idx - 1];
            if (phaseOneCorrectionX < ratio) {
              ratio = phaseOneCorrectionX;
              ixlast = idx;
            }
          }

          if (ixlast == 0) {
            solution->state = 1;
          } else {
            activeSetChangeID = -1;
            globalActiveConstrIdx = ixlast;
            subProblemChanged = true;
            removeConstr(workingset, ixlast);
            solution->lambda[ixlast - 1] = 0.0;
          }
        } else {
          ixlast = workingset->nActiveConstr;
          activeSetChangeID = 0;
          globalActiveConstrIdx = workingset->nActiveConstr;
          subProblemChanged = true;
          removeConstr(workingset, workingset->nActiveConstr);
          solution->lambda[ixlast - 1] = 0.0;
        }

        updateFval = false;
        checkStoppingAndUpdateFval(&activeSetChangeID, f, solution, memspace,
          objective, workingset, qrmanager, options_ObjectiveLimit,
          runTimeOptions_MaxIterations, updateFval);
      }
    } else {
      if (!updateFval) {
        solution->fstar = computeFval_ReuseHx(objective,
          memspace->workspace_double, f, solution->xstar);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

/*
 * Arguments    : boolean_T obj_hasLinear
 *                int obj_nvar
 *                double workspace[28]
 *                const double H[9]
 *                const double f[3]
 *                const double x[4]
 * Return Type  : void
 */
static void linearForm_(boolean_T obj_hasLinear, int obj_nvar, double workspace
  [28], const double H[9], const double f[3], const double x[4])
{
  double c;
  int i;
  int i1;
  int ia;
  int iac;
  int ix;
  int iy;
  ix = 0;
  if (obj_hasLinear) {
    if (0 <= obj_nvar - 1) {
      memcpy(&workspace[0], &f[0], obj_nvar * sizeof(double));
    }

    ix = 1;
  }

  if (obj_nvar != 0) {
    if ((ix != 1) && (0 <= obj_nvar - 1)) {
      memset(&workspace[0], 0, obj_nvar * sizeof(double));
    }

    ix = 0;
    i = obj_nvar * (obj_nvar - 1) + 1;
    for (iac = 1; obj_nvar < 0 ? iac >= i : iac <= i; iac += obj_nvar) {
      c = 0.5 * x[ix];
      iy = 0;
      i1 = (iac + obj_nvar) - 1;
      for (ia = iac; ia <= i1; ia++) {
        workspace[iy] += H[ia - 1] * c;
        iy++;
      }

      ix++;
    }
  }
}

/*
 * Arguments    : const d_struct_T *obj
 *                const double x[28]
 *                int ix0
 * Return Type  : double
 */
static double maxConstraintViolation(const d_struct_T *obj, const double x[28],
  int ix0)
{
  double u1;
  double v;
  int idx;
  int mFixed;
  int mLB;
  int mUB;
  mLB = obj->sizes[3];
  mUB = obj->sizes[4];
  mFixed = obj->sizes[0];
  v = 0.0;
  if (obj->sizes[3] > 0) {
    for (idx = 0; idx < mLB; idx++) {
      u1 = -x[(ix0 + obj->indexLB[idx]) - 2] - obj->lb[obj->indexLB[idx] - 1];
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    for (idx = 0; idx < mUB; idx++) {
      u1 = x[(ix0 + obj->indexUB[idx]) - 2] - obj->ub[obj->indexUB[idx] - 1];
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    for (idx = 0; idx < mFixed; idx++) {
      mLB = obj->indexFixed[idx] - 1;
      u1 = fabs(x[(ix0 + mLB) - 1] - obj->ub[mLB]);
      if ((!(v > u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  return v;
}

/*
 * Arguments    : f_struct_T *obj
 *                int LD_offset
 *                int NColsRemain
 *                double REG_PRIMAL
 * Return Type  : void
 */
static void partialColLDL3_(f_struct_T *obj, int LD_offset, int NColsRemain,
  double REG_PRIMAL)
{
  double a;
  int LD_diagOffset;
  int i;
  int i1;
  int i2;
  int i3;
  int ia;
  int ia0;
  int iac;
  int ix;
  int iy;
  int iy0;
  int j;
  int k;
  int lastColC;
  int m;
  int offsetColK;
  int subRows;
  for (k = 0; k < 48; k++) {
    subRows = (NColsRemain - k) - 1;
    lastColC = k << 3;
    LD_diagOffset = (LD_offset + lastColC) - 1;
    for (ix = 0; ix <= subRows; ix++) {
      obj->workspace_[lastColC + ix] = obj->FMat[LD_diagOffset + ix];
    }

    offsetColK = 7 * k;
    for (ix = 0; ix < NColsRemain; ix++) {
      obj->workspace2_[ix] = obj->workspace_[offsetColK + ix];
    }

    if ((NColsRemain != 0) && (k != 0)) {
      ix = LD_offset + k;
      i = 7 * (k - 1) + 1;
      for (iac = 1; iac <= i; iac += 7) {
        iy = 0;
        i1 = (iac + NColsRemain) - 1;
        for (ia = iac; ia <= i1; ia++) {
          obj->workspace2_[iy] += obj->workspace_[ia - 1] * -obj->FMat[ix - 1];
          iy++;
        }

        ix += 7;
      }
    }

    for (ix = 0; ix < NColsRemain; ix++) {
      obj->workspace_[offsetColK + ix] = obj->workspace2_[ix];
    }

    for (ix = 0; ix <= subRows; ix++) {
      obj->FMat[LD_diagOffset + ix] = obj->workspace_[lastColC + ix];
    }

    if (fabs(obj->FMat[LD_diagOffset]) <= obj->regTol_) {
      obj->FMat[LD_diagOffset] += REG_PRIMAL;
    }

    a = 1.0 / obj->FMat[LD_diagOffset];
    offsetColK = LD_diagOffset + 2;
    i = (LD_diagOffset + subRows) + 1;
    for (ix = offsetColK; ix <= i; ix++) {
      obj->FMat[ix - 1] *= a;
    }
  }

  i = NColsRemain - 1;
  for (j = 48; j <= i; j += 48) {
    subRows = NColsRemain - j;
    if (48 < subRows) {
      LD_diagOffset = 48;
    } else {
      LD_diagOffset = subRows;
    }

    ia0 = j + LD_diagOffset;
    i1 = ia0 - 1;
    for (k = j; k <= i1; k++) {
      m = ia0 - k;
      iy0 = (LD_offset + (k << 3)) - 1;
      for (ix = 0; ix < 48; ix++) {
        obj->workspace2_[ix] = obj->FMat[48];
      }

      offsetColK = k + 1;
      if (m != 0) {
        ix = 0;
        i2 = k + 330;
        for (iac = offsetColK; iac <= i2; iac += 7) {
          iy = iy0;
          i3 = (iac + m) - 1;
          for (ia = iac; ia <= i3; ia++) {
            obj->FMat[iy] += obj->workspace_[ia - 1] * -obj->workspace2_[ix];
            iy++;
          }

          ix++;
        }
      }
    }

    if (ia0 < NColsRemain) {
      m = subRows - LD_diagOffset;
      iac = ((LD_offset + LD_diagOffset) + (j << 3)) - 1;
      for (lastColC = 0; lastColC < 48; lastColC++) {
        offsetColK = (LD_offset + j) + lastColC * 7;
        iy0 = lastColC * 7;
        for (k = 0; k < LD_diagOffset; k++) {
          obj->workspace2_[iy0 + k] = obj->FMat[(offsetColK + k) - 1];
        }
      }

      if ((m != 0) && (LD_diagOffset != 0)) {
        lastColC = iac + 7 * (LD_diagOffset - 1);
        offsetColK = 0;
        for (LD_diagOffset = iac; LD_diagOffset <= lastColC; LD_diagOffset += 7)
        {
          subRows = ia0 - 1;
          offsetColK++;
          i1 = offsetColK + 329;
          for (ix = offsetColK; ix <= i1; ix += 7) {
            ia = subRows;
            i2 = LD_diagOffset + 1;
            i3 = LD_diagOffset + m;
            for (iy = i2; iy <= i3; iy++) {
              ia++;
              obj->FMat[iy - 1] += -obj->workspace2_[ix - 1] * obj->
                workspace_[ia];
            }

            subRows += 7;
          }
        }
      }
    }
  }
}

/*
 * Arguments    : const double H[9]
 *                const double f[3]
 *                struct_T *solution
 *                c_struct_T *memspace
 *                d_struct_T *workingset
 *                e_struct_T *qrmanager
 *                h_struct_T *options
 *                const b_struct_T *runTimeOptions
 *                f_struct_T *cholmanager
 *                g_struct_T *objective
 * Return Type  : void
 */
static void phaseone(const double H[9], const double f[3], struct_T *solution,
                     c_struct_T *memspace, d_struct_T *workingset, e_struct_T
                     *qrmanager, h_struct_T *options, const b_struct_T
                     *runTimeOptions, f_struct_T *cholmanager, g_struct_T
                     *objective)
{
  double d;
  double minLambda;
  int TYPE;
  int activeSetChangeID;
  int b_nVar;
  int exitg1;
  int globalActiveConstrIdx;
  int i;
  int idx;
  int idxEndIneq;
  int idxStartIneq;
  int nVar;
  boolean_T exitg2;
  boolean_T guard1 = false;
  boolean_T subProblemChanged;
  boolean_T updateFval;
  nVar = workingset->nVar;
  solution->xstar[3] = solution->maxConstr + 1.0;
  idxStartIneq = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
  idxEndIneq = workingset->nActiveConstr;
  for (TYPE = idxStartIneq; TYPE <= idxEndIneq; TYPE++) {
    workingset->isActiveConstr[(workingset->isActiveIdx[workingset->Wid[TYPE - 1]
      - 1] + workingset->Wlocalidx[TYPE - 1]) - 2] = false;
  }

  workingset->nWConstr[2] = 0;
  workingset->nWConstr[3] = 0;
  workingset->nWConstr[4] = 0;
  workingset->nActiveConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  setProblemType(workingset, 1);
  objective->grad[0] = 0.0;
  objective->grad[1] = 0.0;
  objective->grad[2] = 0.0;
  objective->grad[3] = 0.0;
  objective->Hx[0] = 0.0;
  objective->Hx[1] = 0.0;
  objective->Hx[2] = 0.0;
  objective->maxVar = 4;
  objective->beta = 0.0;
  objective->rho = 0.0;
  objective->prev_objtype = 3;
  objective->prev_nvar = 3;
  objective->prev_hasLinear = true;
  objective->objtype = 5;
  objective->nvar = 4;
  objective->gammaScalar = 1.0;
  objective->hasLinear = true;
  options->ObjectiveLimit = 1.0E-8;
  options->StepTolerance = 1.4901161193847657E-10;
  solution->fstar = computeFval(objective, memspace->workspace_double, H, f,
    solution->xstar);
  memset(&cholmanager->FMat[0], 0, 49U * sizeof(double));
  cholmanager->ldm = 7;
  cholmanager->ndims = 0;
  cholmanager->info = 0;
  cholmanager->scaleFactor = 100.0;
  cholmanager->ConvexCheck = true;
  cholmanager->regTol_ = 0.0;
  memset(&cholmanager->workspace_[0], 0, 336U * sizeof(double));
  memset(&cholmanager->workspace2_[0], 0, 336U * sizeof(double));
  subProblemChanged = true;
  updateFval = true;
  activeSetChangeID = 0;
  b_nVar = workingset->nVar;
  globalActiveConstrIdx = 0;
  computeGrad_StoreHx(objective, H, f, solution->xstar);
  solution->fstar = computeFval_ReuseHx(objective, memspace->workspace_double, f,
    solution->xstar);
  solution->state = -5;
  for (idxEndIneq = 0; idxEndIneq < 7; idxEndIneq++) {
    solution->lambda[idxEndIneq] = 0.0;
  }

  do {
    exitg1 = 0;
    if (solution->state == -5) {
      guard1 = false;
      if (subProblemChanged) {
        switch (activeSetChangeID) {
         case 1:
          squareQ_appendCol(qrmanager, workingset->ATwset,
                            ((workingset->nActiveConstr - 1) << 2) + 1);
          break;

         case -1:
          deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
          break;

         default:
          factorQR(qrmanager, workingset->ATwset, b_nVar,
                   workingset->nActiveConstr);
          computeQ_(qrmanager, qrmanager->mrows);
          break;
        }

        compute_deltax(H, solution, memspace, qrmanager, cholmanager, objective);
        if (solution->state != -5) {
          exitg1 = 1;
        } else if ((b_xnrm2(b_nVar, solution->searchDir) <
                    1.4901161193847657E-10) || (workingset->nActiveConstr >=
                    b_nVar)) {
          guard1 = true;
        } else {
          feasibleratiotest(solution->xstar, solution->searchDir,
                            workingset->nVar, workingset->lb, workingset->ub,
                            workingset->indexLB, workingset->indexUB,
                            workingset->sizes, workingset->isActiveIdx,
                            workingset->isActiveConstr, workingset->nWConstr,
                            true, &minLambda, &updateFval, &i, &idxStartIneq);
          if (updateFval) {
            switch (i) {
             case 3:
              workingset->nWConstr[2]++;
              workingset->isActiveConstr[(workingset->isActiveIdx[2] +
                idxStartIneq) - 2] = true;
              workingset->nActiveConstr++;
              workingset->Wid[workingset->nActiveConstr - 1] = 3;
              workingset->Wlocalidx[workingset->nActiveConstr - 1] =
                idxStartIneq;

              /* A check that is always false is detected at compile-time. Eliminating code that follows. */
              break;

             case 4:
              addBoundToActiveSetMatrix_(workingset, 4, idxStartIneq);
              break;

             default:
              addBoundToActiveSetMatrix_(workingset, 5, idxStartIneq);
              break;
            }

            activeSetChangeID = 1;
          } else {
            if (objective->objtype == 5) {
              if (b_xnrm2(objective->nvar, solution->searchDir) > 100.0 *
                  (double)objective->nvar * 1.4901161193847656E-8) {
                solution->state = 3;
              } else {
                solution->state = 4;
              }
            }

            subProblemChanged = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }

          if (!(minLambda == 0.0)) {
            idxStartIneq = b_nVar - 1;
            for (idxEndIneq = 0; idxEndIneq <= idxStartIneq; idxEndIneq++) {
              solution->xstar[idxEndIneq] += minLambda * solution->
                searchDir[idxEndIneq];
            }
          }

          computeGrad_StoreHx(objective, H, f, solution->xstar);
          updateFval = true;
          checkStoppingAndUpdateFval(&activeSetChangeID, f, solution, memspace,
            objective, workingset, qrmanager, options->ObjectiveLimit,
            runTimeOptions->MaxIterations, updateFval);
        }
      } else {
        if (0 <= b_nVar - 1) {
          memset(&solution->searchDir[0], 0, b_nVar * sizeof(double));
        }

        guard1 = true;
      }

      if (guard1) {
        compute_lambda(memspace->workspace_double, solution, objective,
                       qrmanager);
        if (solution->state != -7) {
          idxEndIneq = -1;
          minLambda = 0.0 * runTimeOptions->ProbRelTolFactor * 0.0;
          i = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          idxStartIneq = workingset->nActiveConstr;
          for (idx = i; idx <= idxStartIneq; idx++) {
            d = solution->lambda[idx - 1];
            if (d < minLambda) {
              minLambda = d;
              idxEndIneq = idx - 1;
            }
          }

          if (idxEndIneq + 1 == 0) {
            solution->state = 1;
          } else {
            activeSetChangeID = -1;
            globalActiveConstrIdx = idxEndIneq + 1;
            subProblemChanged = true;
            TYPE = workingset->Wid[idxEndIneq] - 1;
            workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
              Wid[idxEndIneq] - 1] + workingset->Wlocalidx[idxEndIneq]) - 2] =
              false;
            workingset->Wid[idxEndIneq] = workingset->Wid
              [workingset->nActiveConstr - 1];
            workingset->Wlocalidx[idxEndIneq] = workingset->Wlocalidx
              [workingset->nActiveConstr - 1];
            i = workingset->nVar;
            for (idx = 0; idx < i; idx++) {
              workingset->ATwset[idx + (idxEndIneq << 2)] = workingset->
                ATwset[idx + ((workingset->nActiveConstr - 1) << 2)];
            }

            workingset->bwset[idxEndIneq] = workingset->bwset
              [workingset->nActiveConstr - 1];
            workingset->nActiveConstr--;
            workingset->nWConstr[TYPE]--;
            solution->lambda[idxEndIneq] = 0.0;
          }
        } else {
          idxEndIneq = workingset->nActiveConstr;
          activeSetChangeID = 0;
          globalActiveConstrIdx = workingset->nActiveConstr;
          subProblemChanged = true;
          idxStartIneq = workingset->Wid[workingset->nActiveConstr - 1] - 1;
          workingset->isActiveConstr[(workingset->isActiveIdx[idxStartIneq] +
            workingset->Wlocalidx[workingset->nActiveConstr - 1]) - 2] = false;
          workingset->nActiveConstr--;
          workingset->nWConstr[idxStartIneq]--;
          solution->lambda[idxEndIneq - 1] = 0.0;
        }

        updateFval = false;
        checkStoppingAndUpdateFval(&activeSetChangeID, f, solution, memspace,
          objective, workingset, qrmanager, options->ObjectiveLimit,
          runTimeOptions->MaxIterations, updateFval);
      }
    } else {
      if (!updateFval) {
        solution->fstar = computeFval_ReuseHx(objective,
          memspace->workspace_double, f, solution->xstar);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (workingset->isActiveConstr[(workingset->isActiveIdx[3] + workingset->
       sizes[3]) - 2]) {
    idx = workingset->sizes[0];
    exitg2 = false;
    while ((!exitg2) && (idx + 1 <= workingset->nActiveConstr)) {
      if ((workingset->Wid[idx] == 4) && (workingset->Wlocalidx[idx] ==
           workingset->sizes[3])) {
        TYPE = workingset->Wid[idx] - 1;
        workingset->isActiveConstr[(workingset->isActiveIdx[workingset->Wid[idx]
          - 1] + workingset->Wlocalidx[idx]) - 2] = false;
        workingset->Wid[idx] = workingset->Wid[workingset->nActiveConstr - 1];
        workingset->Wlocalidx[idx] = workingset->Wlocalidx
          [workingset->nActiveConstr - 1];
        i = workingset->nVar;
        for (idxStartIneq = 0; idxStartIneq < i; idxStartIneq++) {
          workingset->ATwset[idxStartIneq + (idx << 2)] = workingset->
            ATwset[idxStartIneq + ((workingset->nActiveConstr - 1) << 2)];
        }

        workingset->bwset[idx] = workingset->bwset[workingset->nActiveConstr - 1];
        workingset->nActiveConstr--;
        workingset->nWConstr[TYPE]--;
        exitg2 = true;
      } else {
        idx++;
      }
    }
  }

  idxStartIneq = workingset->nActiveConstr - 1;
  while ((idxStartIneq + 1 > workingset->sizes[0]) && (idxStartIneq + 1 > nVar))
  {
    TYPE = workingset->Wid[idxStartIneq] - 1;
    workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
      Wid[idxStartIneq] - 1] + workingset->Wlocalidx[idxStartIneq]) - 2] = false;
    workingset->Wid[idxStartIneq] = workingset->Wid[workingset->nActiveConstr -
      1];
    workingset->Wlocalidx[idxStartIneq] = workingset->Wlocalidx
      [workingset->nActiveConstr - 1];
    i = workingset->nVar;
    for (idx = 0; idx < i; idx++) {
      workingset->ATwset[idx + (idxStartIneq << 2)] = workingset->ATwset[idx +
        ((workingset->nActiveConstr - 1) << 2)];
    }

    workingset->bwset[idxStartIneq] = workingset->bwset
      [workingset->nActiveConstr - 1];
    workingset->nActiveConstr--;
    workingset->nWConstr[TYPE]--;
    idxStartIneq--;
  }

  solution->maxConstr = solution->xstar[3];
  setProblemType(workingset, 3);
  objective->objtype = objective->prev_objtype;
  objective->nvar = objective->prev_nvar;
  objective->hasLinear = objective->prev_hasLinear;
  options->ObjectiveLimit = -1.0E+20;
  options->StepTolerance = 1.0E-8;
}

/*
 * Arguments    : double A[49]
 *                int m
 *                int n
 *                int nfxd
 *                double tau[7]
 * Return Type  : void
 */
static void qrf(double A[49], int m, int n, int nfxd, double tau[7])
{
  double work[7];
  double atmp;
  int i;
  int ii;
  int mmi;
  for (i = 0; i < 7; i++) {
    work[i] = 0.0;
  }

  for (i = 0; i < nfxd; i++) {
    ii = i * 7 + i;
    mmi = m - i;
    if (i + 1 < m) {
      atmp = A[ii];
      tau[i] = xzlarfg(mmi, &atmp, A, ii + 2);
      A[ii] = atmp;
    } else {
      tau[i] = 0.0;
    }

    if (i + 1 < n) {
      atmp = A[ii];
      A[ii] = 1.0;
      xzlarf(mmi, (n - i) - 1, ii + 1, tau[i], A, ii + 8, work);
      A[ii] = atmp;
    }
  }
}

/*
 * Arguments    : d_struct_T *obj
 *                int idx_global
 * Return Type  : void
 */
static void removeConstr(d_struct_T *obj, int idx_global)
{
  int TYPE_tmp;
  int i;
  int idx;
  TYPE_tmp = obj->Wid[idx_global - 1] - 1;
  obj->isActiveConstr[(obj->isActiveIdx[TYPE_tmp] + obj->Wlocalidx[idx_global -
                       1]) - 2] = false;
  obj->Wid[idx_global - 1] = obj->Wid[obj->nActiveConstr - 1];
  obj->Wlocalidx[idx_global - 1] = obj->Wlocalidx[obj->nActiveConstr - 1];
  i = obj->nVar;
  for (idx = 0; idx < i; idx++) {
    obj->ATwset[idx + ((idx_global - 1) << 2)] = obj->ATwset[idx +
      ((obj->nActiveConstr - 1) << 2)];
  }

  obj->bwset[idx_global - 1] = obj->bwset[obj->nActiveConstr - 1];
  obj->nActiveConstr--;
  obj->nWConstr[TYPE_tmp]--;
}

/*
 * Arguments    : d_struct_T *obj
 *                int idx_global
 * Return Type  : void
 */
static void removeEqConstr(d_struct_T *obj, int idx_global)
{
  int TYPE_tmp_tmp;
  int i;
  int idx;
  int totalEq;
  totalEq = (obj->nWConstr[0] + obj->nWConstr[1]) - 1;
  if ((totalEq + 1 != 0) && (idx_global <= totalEq + 1)) {
    if ((obj->nActiveConstr == totalEq + 1) || (idx_global == totalEq + 1)) {
      obj->mEqRemoved++;
      totalEq = obj->Wid[idx_global - 1] - 1;
      obj->isActiveConstr[(obj->isActiveIdx[totalEq] + obj->Wlocalidx[idx_global
                           - 1]) - 2] = false;
      obj->Wid[idx_global - 1] = obj->Wid[obj->nActiveConstr - 1];
      obj->Wlocalidx[idx_global - 1] = obj->Wlocalidx[obj->nActiveConstr - 1];
      i = obj->nVar;
      for (idx = 0; idx < i; idx++) {
        obj->ATwset[idx + ((idx_global - 1) << 2)] = obj->ATwset[idx +
          ((obj->nActiveConstr - 1) << 2)];
      }

      obj->bwset[idx_global - 1] = obj->bwset[obj->nActiveConstr - 1];
      obj->nActiveConstr--;
      obj->nWConstr[totalEq]--;
    } else {
      obj->mEqRemoved++;
      TYPE_tmp_tmp = obj->Wid[idx_global - 1] - 1;
      obj->isActiveConstr[(obj->isActiveIdx[TYPE_tmp_tmp] + obj->
                           Wlocalidx[idx_global - 1]) - 2] = false;
      obj->Wid[idx_global - 1] = obj->Wid[totalEq];
      obj->Wlocalidx[idx_global - 1] = obj->Wlocalidx[totalEq];
      i = obj->nVar;
      for (idx = 0; idx < i; idx++) {
        obj->ATwset[idx + ((idx_global - 1) << 2)] = obj->ATwset[idx + (totalEq <<
          2)];
      }

      obj->bwset[idx_global - 1] = obj->bwset[totalEq];
      obj->Wid[totalEq] = obj->Wid[obj->nActiveConstr - 1];
      obj->Wlocalidx[totalEq] = obj->Wlocalidx[obj->nActiveConstr - 1];
      i = obj->nVar;
      for (idx = 0; idx < i; idx++) {
        obj->ATwset[idx + (totalEq << 2)] = obj->ATwset[idx +
          ((obj->nActiveConstr - 1) << 2)];
      }

      obj->bwset[totalEq] = obj->bwset[obj->nActiveConstr - 1];
      obj->nActiveConstr--;
      obj->nWConstr[TYPE_tmp_tmp]--;
    }
  }
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double y;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = a * sqrt(y * y + 1.0);
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/*
 * Arguments    : d_struct_T *obj
 *                int PROBLEM_TYPE
 * Return Type  : void
 */
static void setProblemType(d_struct_T *obj, int PROBLEM_TYPE)
{
  int i;
  int i1;
  int idx;
  int idxStartIneq;
  int idx_lb;
  switch (PROBLEM_TYPE) {
   case 3:
    obj->nVar = 3;
    obj->mConstr = obj->mConstrOrig;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesNormal[i];
    }

    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxNormal[i];
    }
    break;

   case 1:
    obj->nVar = 4;
    obj->mConstr = obj->mConstrOrig + 1;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesPhaseOne[i];
    }

    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxPhaseOne[i];
    }

    i = obj->sizes[0];
    for (idx = 0; idx < i; idx++) {
      obj->ATwset[(idx << 2) + 3] = 0.0;
    }

    obj->indexLB[obj->sizes[3] - 1] = 4;
    obj->lb[3] = obj->SLACK0;
    idxStartIneq = obj->isActiveIdx[2];
    i = obj->nActiveConstr;
    for (idx = idxStartIneq; idx <= i; idx++) {
      obj->ATwset[((idx - 1) << 2) + 3] = -1.0;
    }
    break;

   case 2:
    obj->nVar = 3;
    obj->mConstr = 6;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesRegularized[i];
    }

    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxRegularized[i];
    }

    if (obj->probType != 4) {
      idx_lb = 3;
      i = obj->sizesNormal[3] + 1;
      i1 = obj->sizesRegularized[3];
      for (idx = i; idx <= i1; idx++) {
        idx_lb++;
        obj->indexLB[idx - 1] = idx_lb;
      }

      idxStartIneq = obj->isActiveIdx[2];
      i = obj->nActiveConstr;
      for (idx = idxStartIneq; idx <= i; idx++) {
        switch (obj->Wid[idx - 1]) {
         case 3:
          i1 = obj->Wlocalidx[idx - 1] + 2;
          for (idx_lb = 4; idx_lb <= i1; idx_lb++) {
            obj->ATwset[((idx - 1) << 2) + 3] = 0.0;
          }

          i1 = (idx - 1) << 2;
          obj->ATwset[(obj->Wlocalidx[idx - 1] + i1) + 2] = -1.0;
          idx_lb = obj->Wlocalidx[idx - 1] + 4;
          if (idx_lb <= 3) {
            memset(&obj->ATwset[(idx_lb + i1) + -1], 0, (((i1 - idx_lb) - i1) +
                    4) * sizeof(double));
          }
          break;
        }
      }
    }
    break;

   default:
    obj->nVar = 4;
    obj->mConstr = 7;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesRegPhaseOne[i];
    }

    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxRegPhaseOne[i];
    }

    i = obj->sizes[0];
    for (idx = 0; idx < i; idx++) {
      obj->ATwset[(idx << 2) + 3] = 0.0;
    }

    obj->indexLB[obj->sizes[3] - 1] = 4;
    obj->lb[3] = obj->SLACK0;
    idxStartIneq = obj->isActiveIdx[2];
    i = obj->nActiveConstr;
    for (idx = idxStartIneq; idx <= i; idx++) {
      obj->ATwset[((idx - 1) << 2) + 3] = -1.0;
    }
    break;
  }

  obj->probType = PROBLEM_TYPE;
}

/*
 * Arguments    : const f_struct_T *obj
 *                double rhs[4]
 * Return Type  : void
 */
static void solve(const f_struct_T *obj, double rhs[4])
{
  double temp;
  int b_i;
  int i;
  int ix;
  int j;
  int jjA;
  int n;
  n = obj->ndims - 2;
  if (obj->ndims != 0) {
    for (j = 0; j <= n + 1; j++) {
      jjA = j + j * 7;
      i = n - j;
      for (b_i = 0; b_i <= i; b_i++) {
        ix = (j + b_i) + 1;
        rhs[ix] -= rhs[j] * obj->FMat[(jjA + b_i) + 1];
      }
    }
  }

  i = obj->ndims;
  for (jjA = 0; jjA < i; jjA++) {
    rhs[jjA] /= obj->FMat[jjA + 7 * jjA];
  }

  n = obj->ndims;
  if (obj->ndims != 0) {
    for (j = n; j >= 1; j--) {
      jjA = (j - 1) * 7;
      temp = rhs[j - 1];
      i = j + 1;
      for (b_i = n; b_i >= i; b_i--) {
        temp -= obj->FMat[(jjA + b_i) - 1] * rhs[b_i - 1];
      }

      rhs[j - 1] = temp;
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double vec[28]
 *                int iv0
 * Return Type  : void
 */
static void squareQ_appendCol(e_struct_T *obj, const double vec[28], int iv0)
{
  double c;
  double s;
  double temp;
  int ia;
  int iac;
  int ix;
  int iy;
  int iyend;
  int n;
  iyend = obj->mrows;
  ix = obj->ncols + 1;
  if (iyend < ix) {
    ix = iyend;
  }

  obj->minRowCol = ix;
  iy = 7 * obj->ncols;
  if (obj->mrows != 0) {
    iyend = iy + obj->mrows;
    if (iy + 1 <= iyend) {
      memset(&obj->QR[iy], 0, (iyend - iy) * sizeof(double));
    }

    iyend = 7 * (obj->mrows - 1) + 1;
    for (iac = 1; iac <= iyend; iac += 7) {
      ix = iv0;
      c = 0.0;
      n = (iac + obj->mrows) - 1;
      for (ia = iac; ia <= n; ia++) {
        c += obj->Q[ia - 1] * vec[ix - 1];
        ix++;
      }

      obj->QR[iy] += c;
      iy++;
    }
  }

  obj->ncols++;
  obj->jpvt[obj->ncols - 1] = obj->ncols;
  for (ix = obj->mrows - 1; ix + 1 > obj->ncols; ix--) {
    temp = obj->QR[ix + 7 * (obj->ncols - 1)];
    xrotg(&obj->QR[(ix + 7 * (obj->ncols - 1)) - 1], &temp, &c, &s);
    obj->QR[ix + 7 * (obj->ncols - 1)] = temp;
    iyend = 7 * (ix - 1);
    n = obj->mrows;
    if (obj->mrows >= 1) {
      iy = iyend + 7;
      for (iac = 0; iac < n; iac++) {
        temp = c * obj->Q[iyend] + s * obj->Q[iy];
        obj->Q[iy] = c * obj->Q[iy] - s * obj->Q[iyend];
        obj->Q[iyend] = temp;
        iy++;
        iyend++;
      }
    }
  }
}

/*
 * Arguments    : int m
 *                int n
 *                int k
 *                const double A[9]
 *                int lda
 *                const double B[49]
 *                int ib0
 *                double C[28]
 * Return Type  : void
 */
static void xgemm(int m, int n, int k, const double A[9], int lda, const double
                  B[49], int ib0, double C[28])
{
  int ar;
  int br;
  int cr;
  int i;
  int i1;
  int i2;
  int ia;
  int ib;
  int ic;
  int lastColC;
  if ((m != 0) && (n != 0)) {
    br = ib0;
    lastColC = 7 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 7) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
      }
    }

    for (cr = 0; cr <= lastColC; cr += 7) {
      ar = -1;
      i = br + k;
      for (ib = br; ib < i; ib++) {
        ia = ar;
        i1 = cr + 1;
        i2 = cr + m;
        for (ic = i1; ic <= i2; ic++) {
          ia++;
          C[ic - 1] += B[ib - 1] * A[ia];
        }

        ar += lda;
      }

      br += 7;
    }
  }
}

/*
 * Arguments    : double A[49]
 *                int m
 *                int n
 *                int jpvt[7]
 *                double tau[7]
 * Return Type  : void
 */
static void xgeqp3(double A[49], int m, int n, int jpvt[7], double tau[7])
{
  double vn1[7];
  double vn2[7];
  double work[7];
  double d;
  double s;
  double temp;
  int b_i;
  int i;
  int ii;
  int ip1;
  int ix;
  int iy;
  int k;
  int minmn_tmp;
  int mmi;
  int nfxd;
  int nmi;
  int pvt;
  if (m < n) {
    minmn_tmp = m;
  } else {
    minmn_tmp = n;
  }

  for (i = 0; i < 7; i++) {
    tau[i] = 0.0;
  }

  if (minmn_tmp < 1) {
    for (pvt = 0; pvt < n; pvt++) {
      jpvt[pvt] = pvt + 1;
    }
  } else {
    nfxd = 0;
    for (pvt = 0; pvt < n; pvt++) {
      if (jpvt[pvt] != 0) {
        nfxd++;
        if (pvt + 1 != nfxd) {
          ix = pvt * 7;
          iy = (nfxd - 1) * 7;
          for (k = 0; k < m; k++) {
            temp = A[ix];
            A[ix] = A[iy];
            A[iy] = temp;
            ix++;
            iy++;
          }

          jpvt[pvt] = jpvt[nfxd - 1];
          jpvt[nfxd - 1] = pvt + 1;
        } else {
          jpvt[pvt] = pvt + 1;
        }
      } else {
        jpvt[pvt] = pvt + 1;
      }
    }

    if (nfxd >= minmn_tmp) {
      nfxd = minmn_tmp;
    }

    qrf(A, m, n, nfxd, tau);
    if (nfxd < minmn_tmp) {
      for (i = 0; i < 7; i++) {
        work[i] = 0.0;
        vn1[i] = 0.0;
        vn2[i] = 0.0;
      }

      b_i = nfxd + 1;
      for (pvt = b_i; pvt <= n; pvt++) {
        d = xnrm2(m - nfxd, A, (nfxd + (pvt - 1) * 7) + 1);
        vn1[pvt - 1] = d;
        vn2[pvt - 1] = d;
      }

      b_i = nfxd + 1;
      for (i = b_i; i <= minmn_tmp; i++) {
        ip1 = i + 1;
        iy = (i - 1) * 7;
        ii = (iy + i) - 1;
        nmi = (n - i) + 1;
        mmi = m - i;
        if (nmi < 1) {
          nfxd = -2;
        } else {
          nfxd = -1;
          if (nmi > 1) {
            ix = i;
            temp = fabs(vn1[i - 1]);
            for (k = 2; k <= nmi; k++) {
              ix++;
              s = fabs(vn1[ix - 1]);
              if (s > temp) {
                nfxd = k - 2;
                temp = s;
              }
            }
          }
        }

        pvt = i + nfxd;
        if (pvt + 1 != i) {
          ix = pvt * 7;
          for (k = 0; k < m; k++) {
            temp = A[ix];
            A[ix] = A[iy];
            A[iy] = temp;
            ix++;
            iy++;
          }

          nfxd = jpvt[pvt];
          jpvt[pvt] = jpvt[i - 1];
          jpvt[i - 1] = nfxd;
          vn1[pvt] = vn1[i - 1];
          vn2[pvt] = vn2[i - 1];
        }

        if (i < m) {
          temp = A[ii];
          d = xzlarfg(mmi + 1, &temp, A, ii + 2);
          tau[i - 1] = d;
          A[ii] = temp;
        } else {
          d = 0.0;
          tau[i - 1] = 0.0;
        }

        if (i < n) {
          temp = A[ii];
          A[ii] = 1.0;
          xzlarf(mmi + 1, nmi - 1, ii + 1, d, A, ii + 8, work);
          A[ii] = temp;
        }

        for (pvt = ip1; pvt <= n; pvt++) {
          nfxd = i + (pvt - 1) * 7;
          d = vn1[pvt - 1];
          if (d != 0.0) {
            temp = fabs(A[nfxd - 1]) / d;
            temp = 1.0 - temp * temp;
            if (temp < 0.0) {
              temp = 0.0;
            }

            s = d / vn2[pvt - 1];
            s = temp * (s * s);
            if (s <= 1.4901161193847656E-8) {
              if (i < m) {
                d = xnrm2(mmi, A, nfxd + 1);
                vn1[pvt - 1] = d;
                vn2[pvt - 1] = d;
              } else {
                vn1[pvt - 1] = 0.0;
                vn2[pvt - 1] = 0.0;
              }
            } else {
              vn1[pvt - 1] = d * sqrt(temp);
            }
          }
        }
      }
    }
  }
}

/*
 * Arguments    : int n
 *                const double x[49]
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const double x[49], int ix0)
{
  double absxk;
  double scale;
  double t;
  double y;
  int k;
  int kend;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : double *a
 *                double *b
 *                double *c
 *                double *s
 * Return Type  : void
 */
static void xrotg(double *a, double *b, double *c, double *s)
{
  double absa;
  double absb;
  double ads;
  double bds;
  double roe;
  double scale;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    ads = absa / scale;
    bds = absb / scale;
    scale *= sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }

    *c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = scale;
  }
}

/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                double tau
 *                double C[49]
 *                int ic0
 *                double work[7]
 * Return Type  : void
 */
static void xzlarf(int m, int n, int iv0, double tau, double C[49], int ic0,
                   double work[7])
{
  double c;
  int b_i;
  int exitg1;
  int i;
  int ia;
  int iac;
  int ix;
  int jy;
  int lastc;
  int lastv;
  boolean_T exitg2;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * 7;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      if (0 <= lastc) {
        memset(&work[0], 0, (lastc + 1) * sizeof(double));
      }

      i = 0;
      b_i = ic0 + 7 * lastc;
      for (iac = ic0; iac <= b_i; iac += 7) {
        ix = iv0;
        c = 0.0;
        jy = (iac + lastv) - 1;
        for (ia = iac; ia <= jy; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0;
      jy = 0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += 7;
      }
    }
  }
}

/*
 * Arguments    : int n
 *                double *alpha1
 *                double x[49]
 *                int ix0
 * Return Type  : double
 */
static double xzlarfg(int n, double *alpha1, double x[49], int ix0)
{
  double beta1;
  double tau;
  double xnorm;
  int i;
  int k;
  int knt;
  tau = 0.0;
  if (n > 0) {
    xnorm = xnrm2(n - 1, x, ix0);
    if (xnorm != 0.0) {
      beta1 = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        i = (ix0 + n) - 2;
        do {
          knt++;
          for (k = ix0; k <= i; k++) {
            x[k - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(beta1) >= 1.0020841800044864E-292));

        beta1 = rt_hypotd_snf(*alpha1, xnrm2(n - 1, x, ix0));
        if (*alpha1 >= 0.0) {
          beta1 = -beta1;
        }

        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        for (k = ix0; k <= i; k++) {
          x[k - 1] *= xnorm;
        }

        for (k = 0; k <= knt; k++) {
          beta1 *= 1.0020841800044864E-292;
        }

        *alpha1 = beta1;
      } else {
        tau = (beta1 - *alpha1) / beta1;
        xnorm = 1.0 / (*alpha1 - beta1);
        i = (ix0 + n) - 2;
        for (k = ix0; k <= i; k++) {
          x[k - 1] *= xnorm;
        }

        *alpha1 = beta1;
      }
    }
  }

  return tau;
}

/*
 * H = [1,-1,1
 *      -1,2,-2
 *      1,-2,4];
 *  f = [2;-3;1];
 *  lb = zeros(3,1);
 *  ub = ones(size(lb));
 *  Aeq = ones(1,3);
 *  beq = 1/2;
 * Arguments    : const double H[9]
 *                const double f[3]
 *                const double lb[3]
 *                const double ub[3]
 *                double x[3]
 *                double *fval
 * Return Type  : void
 */
void test_quadp(const double H[9], const double f[3], const double lb[3], const
                double ub[3], double x[3], double *fval)
{
  b_struct_T expl_temp;
  c_struct_T memspace;
  d_struct_T WorkingSet;
  e_struct_T QRManager;
  f_struct_T CholRegManager;
  g_struct_T QPObjective;
  struct_T solution;
  double H_infnrm;
  double colSum;
  double f_infnrm;
  int i;
  int i2;
  int i3;
  int idxFillStart;
  int idx_local;
  int mFixed;
  int mLB;
  int mUB;
  signed char b_obj_tmp[5];
  signed char obj_tmp[5];
  signed char indexFixed[3];
  signed char indexLB[3];
  signed char indexUB[3];
  signed char b_i;
  signed char i1;
  boolean_T guard1 = false;
  if (!isInitialized_test_quadp) {
    test_quadp_initialize();
  }

  solution.fstar = 0.0;
  solution.firstorderopt = 0.0;
  for (i = 0; i < 7; i++) {
    solution.lambda[i] = 0.0;
  }

  solution.state = 0;
  solution.maxConstr = 0.0;
  solution.iterations = 0;
  solution.searchDir[0] = 0.0;
  solution.searchDir[1] = 0.0;
  solution.searchDir[2] = 0.0;
  solution.searchDir[3] = 0.0;
  mLB = 0;
  mUB = 0;
  mFixed = 0;
  solution.xstar[0] = 0.0;
  guard1 = false;
  if ((!rtIsInf(lb[0])) && (!rtIsNaN(lb[0]))) {
    if (fabs(lb[0] - ub[0]) < 1.0E-8) {
      mFixed = 1;
      indexFixed[0] = 1;
    } else {
      mLB = 1;
      indexLB[0] = 1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 && ((!rtIsInf(ub[0])) && (!rtIsNaN(ub[0])))) {
    mUB = 1;
    indexUB[0] = 1;
  }

  solution.xstar[1] = 0.0;
  guard1 = false;
  if ((!rtIsInf(lb[1])) && (!rtIsNaN(lb[1]))) {
    if (fabs(lb[1] - ub[1]) < 1.0E-8) {
      mFixed++;
      indexFixed[mFixed - 1] = 2;
    } else {
      mLB++;
      indexLB[mLB - 1] = 2;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 && ((!rtIsInf(ub[1])) && (!rtIsNaN(ub[1])))) {
    mUB++;
    indexUB[mUB - 1] = 2;
  }

  solution.xstar[2] = 0.0;
  guard1 = false;
  if ((!rtIsInf(lb[2])) && (!rtIsNaN(lb[2]))) {
    if (fabs(lb[2] - ub[2]) < 1.0E-8) {
      mFixed++;
      indexFixed[mFixed - 1] = 3;
    } else {
      mLB++;
      indexLB[mLB - 1] = 3;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1 && ((!rtIsInf(ub[2])) && (!rtIsNaN(ub[2])))) {
    mUB++;
    indexUB[mUB - 1] = 3;
  }

  i = (mLB + mUB) + mFixed;
  WorkingSet.mConstr = i;
  WorkingSet.mConstrOrig = i;
  WorkingSet.mConstrMax = 7;
  WorkingSet.nVar = 3;
  WorkingSet.nVarOrig = 3;
  WorkingSet.nVarMax = 4;
  WorkingSet.ldA = 4;
  WorkingSet.mEqRemoved = 0;
  WorkingSet.nActiveConstr = 0;
  obj_tmp[0] = (signed char)mFixed;
  obj_tmp[1] = 0;
  obj_tmp[2] = 0;
  obj_tmp[3] = (signed char)mLB;
  obj_tmp[4] = (signed char)mUB;
  b_obj_tmp[0] = (signed char)mFixed;
  b_obj_tmp[1] = 0;
  b_obj_tmp[2] = 0;
  b_obj_tmp[3] = (signed char)(mLB + 1);
  b_obj_tmp[4] = (signed char)mUB;
  WorkingSet.isActiveIdx[0] = 1;
  WorkingSet.isActiveIdx[1] = mFixed;
  WorkingSet.isActiveIdx[2] = 0;
  WorkingSet.isActiveIdx[3] = 0;
  WorkingSet.isActiveIdx[4] = mLB;
  WorkingSet.isActiveIdx[5] = mUB;
  for (i = 0; i < 5; i++) {
    b_i = obj_tmp[i];
    WorkingSet.sizes[i] = b_i;
    WorkingSet.sizesNormal[i] = b_i;
    i1 = b_obj_tmp[i];
    WorkingSet.sizesPhaseOne[i] = i1;
    WorkingSet.sizesRegularized[i] = b_i;
    WorkingSet.sizesRegPhaseOne[i] = i1;
    WorkingSet.isActiveIdx[i + 1] += WorkingSet.isActiveIdx[i];
  }

  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdxNormal[i] = WorkingSet.isActiveIdx[i];
  }

  WorkingSet.isActiveIdxPhaseOne[0] = 1;
  WorkingSet.isActiveIdxPhaseOne[1] = mFixed;
  WorkingSet.isActiveIdxPhaseOne[2] = 0;
  WorkingSet.isActiveIdxPhaseOne[3] = 0;
  WorkingSet.isActiveIdxPhaseOne[4] = mLB + 1;
  WorkingSet.isActiveIdxPhaseOne[5] = mUB;
  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxPhaseOne[i + 1] += WorkingSet.isActiveIdxPhaseOne[i];
  }

  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdxRegularized[i] = WorkingSet.isActiveIdx[i];
    WorkingSet.isActiveIdxRegPhaseOne[i] = WorkingSet.isActiveIdxPhaseOne[i];
  }

  for (i = 0; i < 5; i++) {
    WorkingSet.nWConstr[i] = 0;
  }

  WorkingSet.probType = 3;
  WorkingSet.SLACK0 = 1.0E-5;
  WorkingSet.lb[0] = -lb[0];
  WorkingSet.ub[0] = ub[0];
  WorkingSet.lb[1] = -lb[1];
  WorkingSet.ub[1] = ub[1];
  WorkingSet.lb[2] = -lb[2];
  WorkingSet.ub[2] = ub[2];
  for (i = 0; i < mLB; i++) {
    WorkingSet.indexLB[i] = indexLB[i];
  }

  for (i = 0; i < mUB; i++) {
    WorkingSet.indexUB[i] = indexUB[i];
  }

  for (i = 0; i < mFixed; i++) {
    WorkingSet.indexFixed[i] = indexFixed[i];
  }

  setProblemType(&WorkingSet, 3);
  idxFillStart = WorkingSet.isActiveIdx[2];
  for (i = idxFillStart; i < 8; i++) {
    WorkingSet.isActiveConstr[i - 1] = false;
  }

  WorkingSet.nWConstr[0] = WorkingSet.sizes[0];
  WorkingSet.nWConstr[1] = 0;
  WorkingSet.nWConstr[2] = 0;
  WorkingSet.nWConstr[3] = 0;
  WorkingSet.nWConstr[4] = 0;
  WorkingSet.nActiveConstr = WorkingSet.nWConstr[0];
  idxFillStart = WorkingSet.sizes[0];
  for (idx_local = 0; idx_local < idxFillStart; idx_local++) {
    WorkingSet.Wid[idx_local] = 1;
    WorkingSet.Wlocalidx[idx_local] = idx_local + 1;
    WorkingSet.isActiveConstr[idx_local] = true;
    i = WorkingSet.indexFixed[idx_local];
    if (0 <= i - 2) {
      memset(&WorkingSet.ATwset[idx_local * 4], 0, (i + -1) * sizeof(double));
    }

    i = idx_local << 2;
    WorkingSet.ATwset[(WorkingSet.indexFixed[idx_local] + i) - 1] = 1.0;
    i2 = WorkingSet.indexFixed[idx_local] + 1;
    i3 = WorkingSet.nVar;
    if (i2 <= i3) {
      memset(&WorkingSet.ATwset[(i2 + i) + -1], 0, ((((i3 + i) - i2) - i) + 1) *
             sizeof(double));
    }

    WorkingSet.bwset[idx_local] = WorkingSet.ub[WorkingSet.indexFixed[idx_local]
      - 1];
  }

  WorkingSet.SLACK0 = 0.0;
  H_infnrm = 0.0;
  f_infnrm = 0.0;
  for (i = 0; i < 3; i++) {
    colSum = (fabs(H[3 * i]) + fabs(H[3 * i + 1])) + fabs(H[3 * i + 2]);
    if ((!(H_infnrm > colSum)) && (!rtIsNaN(colSum))) {
      H_infnrm = colSum;
    }

    colSum = fabs(f[i]);
    if ((!(f_infnrm > colSum)) && (!rtIsNaN(colSum))) {
      f_infnrm = colSum;
    }
  }

  if ((1.0 > f_infnrm) || rtIsNaN(f_infnrm)) {
    colSum = 1.0;
  } else {
    colSum = f_infnrm;
  }

  if ((!(colSum > H_infnrm)) && (!rtIsNaN(H_infnrm))) {
    colSum = H_infnrm;
  }

  if ((1.0 > colSum) || rtIsNaN(colSum)) {
    colSum = 1.0;
  }

  expl_temp.RemainFeasible = false;
  expl_temp.ProbRelTolFactor = colSum;
  expl_temp.ConstrRelTolFactor = 1.0;
  expl_temp.MaxIterations = 10 * (((mFixed + mLB) + mUB) + 3);
  driver(H, f, &solution, &memspace, &WorkingSet, expl_temp, &QRManager,
         &CholRegManager, &QPObjective);
  x[0] = solution.xstar[0];
  x[1] = solution.xstar[1];
  x[2] = solution.xstar[2];
  if (solution.state > 0) {
    *fval = solution.fstar;
  } else {
    *fval = computeFval(&QPObjective, memspace.workspace_double, H, f,
                        solution.xstar);
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void test_quadp_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_test_quadp = true;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void test_quadp_terminate(void)
{
  /* (no terminate code required) */
  isInitialized_test_quadp = false;
}

/*
 * File trailer for test_quadp.c
 *
 * [EOF]
 */
