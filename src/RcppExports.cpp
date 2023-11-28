// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// elasticNetFit
Rcpp::List elasticNetFit(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> coefInit, double lambda1, double lambda2, int maxIter, double tol, bool verbose, int verboseNumIter, int verboseNumCoef);
RcppExport SEXP _HDMAADMM_elasticNetFit(SEXP XSEXP, SEXP ySEXP, SEXP coefInitSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP maxIterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP verboseNumIterSEXP, SEXP verboseNumCoefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type coefInit(coefInitSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumIter(verboseNumIterSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumCoef(verboseNumCoefSEXP);
    rcpp_result_gen = Rcpp::wrap(elasticNetFit(X, y, coefInit, lambda1, lambda2, maxIter, tol, verbose, verboseNumIter, verboseNumCoef));
    return rcpp_result_gen;
END_RCPP
}
// fMatProd
Eigen::MatrixXd fMatProd(SEXP X, SEXP Y, bool is_X_symmetric);
RcppExport SEXP _HDMAADMM_fMatProd(SEXP XSEXP, SEXP YSEXP, SEXP is_X_symmetricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type is_X_symmetric(is_X_symmetricSEXP);
    rcpp_result_gen = Rcpp::wrap(fMatProd(X, Y, is_X_symmetric));
    return rcpp_result_gen;
END_RCPP
}
// fMatTransProd
Eigen::MatrixXd fMatTransProd(SEXP X, SEXP Y, bool is_X_symmetric);
RcppExport SEXP _HDMAADMM_fMatTransProd(SEXP XSEXP, SEXP YSEXP, SEXP is_X_symmetricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type is_X_symmetric(is_X_symmetricSEXP);
    rcpp_result_gen = Rcpp::wrap(fMatTransProd(X, Y, is_X_symmetric));
    return rcpp_result_gen;
END_RCPP
}
// fMatInv
Eigen::MatrixXd fMatInv(SEXP X, bool is_sym_pd);
RcppExport SEXP _HDMAADMM_fMatInv(SEXP XSEXP, SEXP is_sym_pdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sym_pd(is_sym_pdSEXP);
    rcpp_result_gen = Rcpp::wrap(fMatInv(X, is_sym_pd));
    return rcpp_result_gen;
END_RCPP
}
// fMatChol
Eigen::MatrixXd fMatChol(SEXP X);
RcppExport SEXP _HDMAADMM_fMatChol(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fMatChol(X));
    return rcpp_result_gen;
END_RCPP
}
// singleModalityAdmmFit
Rcpp::List singleModalityAdmmFit(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::MatrixXd> M1, Eigen::Map<Eigen::MatrixXd> alphaInit, Eigen::Map<Eigen::MatrixXd> betaInit, Eigen::Map<Eigen::MatrixXd> gammaInit, double rho, double lambda1a, double lambda1b, double lambda1g, double lambda2a, double lambda2b, int penaltyType, Rcpp::List penaltyParameters, Eigen::Map<Eigen::MatrixXd> XtX, Eigen::Map<Eigen::MatrixXd> XtXInv, Eigen::Map<Eigen::MatrixXd> XtXPlusRhoInv, Eigen::Map<Eigen::MatrixXd> XtM1, Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv, Eigen::Map<Eigen::MatrixXd> M1tY, Eigen::Map<Eigen::MatrixXd> XtY, int maxIter, double tol, bool verbose, int verboseNumIter, int verboseNumAlpha, int verboseNumBeta, int verboseNumGamma);
RcppExport SEXP _HDMAADMM_singleModalityAdmmFit(SEXP XSEXP, SEXP YSEXP, SEXP M1SEXP, SEXP alphaInitSEXP, SEXP betaInitSEXP, SEXP gammaInitSEXP, SEXP rhoSEXP, SEXP lambda1aSEXP, SEXP lambda1bSEXP, SEXP lambda1gSEXP, SEXP lambda2aSEXP, SEXP lambda2bSEXP, SEXP penaltyTypeSEXP, SEXP penaltyParametersSEXP, SEXP XtXSEXP, SEXP XtXInvSEXP, SEXP XtXPlusRhoInvSEXP, SEXP XtM1SEXP, SEXP M1tM1PlusRhoInvSEXP, SEXP M1tYSEXP, SEXP XtYSEXP, SEXP maxIterSEXP, SEXP tolSEXP, SEXP verboseSEXP, SEXP verboseNumIterSEXP, SEXP verboseNumAlphaSEXP, SEXP verboseNumBetaSEXP, SEXP verboseNumGammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type alphaInit(alphaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1a(lambda1aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1b(lambda1bSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1g(lambda1gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2a(lambda2aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2b(lambda2bSEXP);
    Rcpp::traits::input_parameter< int >::type penaltyType(penaltyTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type penaltyParameters(penaltyParametersSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtXInv(XtXInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtXPlusRhoInv(XtXPlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtM1(XtM1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tM1PlusRhoInv(M1tM1PlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tY(M1tYSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtY(XtYSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumIter(verboseNumIterSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumAlpha(verboseNumAlphaSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumBeta(verboseNumBetaSEXP);
    Rcpp::traits::input_parameter< int >::type verboseNumGamma(verboseNumGammaSEXP);
    rcpp_result_gen = Rcpp::wrap(singleModalityAdmmFit(X, Y, M1, alphaInit, betaInit, gammaInit, rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b, penaltyType, penaltyParameters, XtX, XtXInv, XtXPlusRhoInv, XtM1, M1tM1PlusRhoInv, M1tY, XtY, maxIter, tol, verbose, verboseNumIter, verboseNumAlpha, verboseNumBeta, verboseNumGamma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HDMAADMM_elasticNetFit", (DL_FUNC) &_HDMAADMM_elasticNetFit, 10},
    {"_HDMAADMM_fMatProd", (DL_FUNC) &_HDMAADMM_fMatProd, 3},
    {"_HDMAADMM_fMatTransProd", (DL_FUNC) &_HDMAADMM_fMatTransProd, 3},
    {"_HDMAADMM_fMatInv", (DL_FUNC) &_HDMAADMM_fMatInv, 2},
    {"_HDMAADMM_fMatChol", (DL_FUNC) &_HDMAADMM_fMatChol, 1},
    {"_HDMAADMM_singleModalityAdmmFit", (DL_FUNC) &_HDMAADMM_singleModalityAdmmFit, 28},
    {NULL, NULL, 0}
};

RcppExport void R_init_HDMAADMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
