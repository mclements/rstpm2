#include <R.h>
#include <Rinternals.h>
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

static const R_CallMethodDef CallEntries[] = {
    {"model_output",       (DL_FUNC) &model_output,       1},
    {"test_cox_tvc2",      (DL_FUNC) &test_cox_tvc2,      1},
    {"test_cox_tvc2_grad", (DL_FUNC) &test_cox_tvc2_grad, 1},
    {"fitCureModel", (DL_FUNC) &fitCureModel, 6},
    {"aft_model_output", (DL_FUNC) &aft_model_output, 1},
    {NULL, NULL, 0}
};

void R_init_rstpm2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
