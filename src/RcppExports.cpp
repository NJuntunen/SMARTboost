// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// loopfeatures_cpp
NumericMatrix loopfeatures_cpp(NumericVector r, NumericVector h, NumericMatrix G0, NumericMatrix x, IntegerVector ifit, IntegerVector infeatures, NumericMatrix mugrid, LogicalVector dichotomous, NumericVector taugrid, List param, double var_epsilon);
RcppExport SEXP _SMARTboost_loopfeatures_cpp(SEXP rSEXP, SEXP hSEXP, SEXP G0SEXP, SEXP xSEXP, SEXP ifitSEXP, SEXP infeaturesSEXP, SEXP mugridSEXP, SEXP dichotomousSEXP, SEXP taugridSEXP, SEXP paramSEXP, SEXP var_epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G0(G0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ifit(ifitSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type infeatures(infeaturesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mugrid(mugridSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type dichotomous(dichotomousSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type taugrid(taugridSEXP);
    Rcpp::traits::input_parameter< List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< double >::type var_epsilon(var_epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(loopfeatures_cpp(r, h, G0, x, ifit, infeatures, mugrid, dichotomous, taugrid, param, var_epsilon));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMARTboost_loopfeatures_cpp", (DL_FUNC) &_SMARTboost_loopfeatures_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMARTboost(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
