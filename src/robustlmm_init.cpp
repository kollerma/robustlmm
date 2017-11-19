#include "globals.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
  
  CALLDEF(_rcpp_module_boot_psi_function_module, 0),
  {NULL, NULL, 0}
  
};

extern "C" void R_unload_robustlmm(DllInfo *info) {  // #nocov start
  // Release resources
} 						// #nocov end

extern "C" void R_init_robustlmm(DllInfo* dllinfo) {
  R_registerRoutines(dllinfo, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dllinfo, FALSE);	// set up symbol symbol lookup (cf R 3.4.0)
}
