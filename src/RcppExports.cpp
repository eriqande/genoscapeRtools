// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mat_cumul_cpp
NumericMatrix mat_cumul_cpp(NumericMatrix x, int dim);
RcppExport SEXP _genoscapeRtools_mat_cumul_cpp(SEXP xSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_cumul_cpp(x, dim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_genoscapeRtools_mat_cumul_cpp", (DL_FUNC) &_genoscapeRtools_mat_cumul_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_genoscapeRtools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
