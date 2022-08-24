#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP model_output(SEXP);
extern SEXP test_cox_tvc2(SEXP);
extern SEXP test_cox_tvc2_grad(SEXP);
extern SEXP fitCureModel(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP aft_model_output(SEXP);
extern SEXP vunirootRcpp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP multistate_ddt(SEXP,SEXP,SEXP,SEXP);
extern SEXP runMarkovODE(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP runExpM(SEXP,SEXP,SEXP);
extern SEXP plugin_P_by(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP plugin_P_L_by(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP test_read_gsm(SEXP);

/* .Fortran calls -- thanks to Gordon Smyth */
extern void F77_NAME(gausq2)(void *, void *, void *, void *, void *);
extern void F77_NAME(laplace_pvf_e1)(void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"model_output",       (DL_FUNC) &model_output,       1},
    {"test_cox_tvc2",      (DL_FUNC) &test_cox_tvc2,      1},
    {"test_cox_tvc2_grad", (DL_FUNC) &test_cox_tvc2_grad, 1},
    {"fitCureModel",       (DL_FUNC) &fitCureModel,       6},
    {"aft_model_output",   (DL_FUNC) &aft_model_output,   1},
    {"vunirootRcpp",       (DL_FUNC) &vunirootRcpp,       7},
    {"multistate_ddt",     (DL_FUNC) &multistate_ddt,     4},
    {"runMarkovODE",       (DL_FUNC) &runMarkovODE,       9},
    {"runExpM",            (DL_FUNC) &runExpM,            3},
    {"plugin_P_by",        (DL_FUNC) &plugin_P_by,        7},
    {"plugin_P_L_by",      (DL_FUNC) &plugin_P_L_by,     10},
    {"test_read_gsm",      (DL_FUNC) &test_read_gsm,      1},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"gausq2", (DL_FUNC) &F77_NAME(gausq2), 5},
    {"laplace_pvf_e1", (DL_FUNC) &F77_NAME(laplace_pvf_e1), 6},
    {NULL, NULL, 0}
};

void R_init_rstpm2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
