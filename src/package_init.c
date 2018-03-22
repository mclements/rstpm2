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

/* .Fortran calls -- thanks for Gordon Smyth */
extern void F77_NAME(gausq2)(void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"model_output",       (DL_FUNC) &model_output,       1},
    {"test_cox_tvc2",      (DL_FUNC) &test_cox_tvc2,      1},
    {"test_cox_tvc2_grad", (DL_FUNC) &test_cox_tvc2_grad, 1},
    {"fitCureModel", (DL_FUNC) &fitCureModel, 6},
    {"aft_model_output", (DL_FUNC) &aft_model_output, 1},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"gausq2", (DL_FUNC) &F77_NAME(gausq2), 5},
    {NULL, NULL, 0}
};

void R_init_rstpm2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
