// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getRegionCPP
IntegerVector getRegionCPP(NumericVector x, NumericVector y);
RcppExport SEXP _rgm_getRegionCPP(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(getRegionCPP(x, y));
    return rcpp_result_gen;
END_RCPP
}
// checkRegionCPP
int checkRegionCPP(int x, int y);
RcppExport SEXP _rgm_checkRegionCPP(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(checkRegionCPP(x, y));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _rgm_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rgm_getRegionCPP", (DL_FUNC) &_rgm_getRegionCPP, 2},
    {"_rgm_checkRegionCPP", (DL_FUNC) &_rgm_checkRegionCPP, 2},
    {"_rgm_timesTwo", (DL_FUNC) &_rgm_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rgm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
