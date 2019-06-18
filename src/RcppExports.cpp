// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getFastqCount
long int getFastqCount();
RcppExport SEXP _nanopoRe_getFastqCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getFastqCount());
    return rcpp_result_gen;
END_RCPP
}
// getFastqBases
long long int getFastqBases();
RcppExport SEXP _nanopoRe_getFastqBases() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getFastqBases());
    return rcpp_result_gen;
END_RCPP
}
// getMalformedFastqHeaderCount
int getMalformedFastqHeaderCount();
RcppExport SEXP _nanopoRe_getMalformedFastqHeaderCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getMalformedFastqHeaderCount());
    return rcpp_result_gen;
END_RCPP
}
// getFastqPlusErrorCount
int getFastqPlusErrorCount();
RcppExport SEXP _nanopoRe_getFastqPlusErrorCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getFastqPlusErrorCount());
    return rcpp_result_gen;
END_RCPP
}
// getZeroLengthSequenceCount
int getZeroLengthSequenceCount();
RcppExport SEXP _nanopoRe_getZeroLengthSequenceCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getZeroLengthSequenceCount());
    return rcpp_result_gen;
END_RCPP
}
// getSequenceQualityMismatchCount
int getSequenceQualityMismatchCount();
RcppExport SEXP _nanopoRe_getSequenceQualityMismatchCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getSequenceQualityMismatchCount());
    return rcpp_result_gen;
END_RCPP
}
// getSkippedLineCount
int getSkippedLineCount();
RcppExport SEXP _nanopoRe_getSkippedLineCount() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getSkippedLineCount());
    return rcpp_result_gen;
END_RCPP
}
// fastqValidator
LogicalVector fastqValidator(std::string fastq);
RcppExport SEXP _nanopoRe_fastqValidator(SEXP fastqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastq(fastqSEXP);
    rcpp_result_gen = Rcpp::wrap(fastqValidator(fastq));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _nanopoRe_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nanopoRe_getFastqCount", (DL_FUNC) &_nanopoRe_getFastqCount, 0},
    {"_nanopoRe_getFastqBases", (DL_FUNC) &_nanopoRe_getFastqBases, 0},
    {"_nanopoRe_getMalformedFastqHeaderCount", (DL_FUNC) &_nanopoRe_getMalformedFastqHeaderCount, 0},
    {"_nanopoRe_getFastqPlusErrorCount", (DL_FUNC) &_nanopoRe_getFastqPlusErrorCount, 0},
    {"_nanopoRe_getZeroLengthSequenceCount", (DL_FUNC) &_nanopoRe_getZeroLengthSequenceCount, 0},
    {"_nanopoRe_getSequenceQualityMismatchCount", (DL_FUNC) &_nanopoRe_getSequenceQualityMismatchCount, 0},
    {"_nanopoRe_getSkippedLineCount", (DL_FUNC) &_nanopoRe_getSkippedLineCount, 0},
    {"_nanopoRe_fastqValidator", (DL_FUNC) &_nanopoRe_fastqValidator, 1},
    {"_nanopoRe_timesTwo", (DL_FUNC) &_nanopoRe_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_nanopoRe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
