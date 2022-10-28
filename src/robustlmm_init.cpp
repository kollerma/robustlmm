#include "globals.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

  CALLDEF(_rcpp_module_boot_psi_function_module, 0),
  CALLDEF(_rcpp_module_boot_rlmerMatrixUtils_module, 0),
  {NULL, NULL, 0}

};

cholmod_common c;

extern "C" void R_unload_robustlmm(DllInfo *info) {
    M_cholmod_finish(&c);
}

extern "C" void R_init_robustlmm(DllInfo* dllinfo) {
  R_registerRoutines(dllinfo, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dllinfo, FALSE);	// set up symbol symbol lookup (cf R 3.4.0)

  M_R_cholmod_start(&c);
}
