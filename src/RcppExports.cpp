// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// estimateElasticNet
Rcpp::List estimateElasticNet(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::MatrixXd> M1, Eigen::Map<Eigen::MatrixXd> alpha, Eigen::Map<Eigen::MatrixXd> beta, Eigen::Map<Eigen::MatrixXd> gamma, Eigen::Map<Eigen::MatrixXd> tauAlpha, Eigen::Map<Eigen::MatrixXd> tauBeta, double rho, double lambda1a, double lambda1b, double lambda1g, double lambda2a, double lambda2b, Eigen::Map<Eigen::MatrixXd> XtX, Eigen::Map<Eigen::MatrixXd> XtXPlusRhoInv, Eigen::Map<Eigen::MatrixXd> XtM1, Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv, Eigen::Map<Eigen::MatrixXd> M1tY, Eigen::Map<Eigen::MatrixXd> XtY);
RcppExport SEXP _HDMAADMM_estimateElasticNet(SEXP XSEXP, SEXP YSEXP, SEXP M1SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauAlphaSEXP, SEXP tauBetaSEXP, SEXP rhoSEXP, SEXP lambda1aSEXP, SEXP lambda1bSEXP, SEXP lambda1gSEXP, SEXP lambda2aSEXP, SEXP lambda2bSEXP, SEXP XtXSEXP, SEXP XtXPlusRhoInvSEXP, SEXP XtM1SEXP, SEXP M1tM1PlusRhoInvSEXP, SEXP M1tYSEXP, SEXP XtYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauAlpha(tauAlphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauBeta(tauBetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1a(lambda1aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1b(lambda1bSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1g(lambda1gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2a(lambda2aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2b(lambda2bSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtXPlusRhoInv(XtXPlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtM1(XtM1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tM1PlusRhoInv(M1tM1PlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tY(M1tYSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtY(XtYSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateElasticNet(X, Y, M1, alpha, beta, gamma, tauAlpha, tauBeta, rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b, XtX, XtXPlusRhoInv, XtM1, M1tM1PlusRhoInv, M1tY, XtY));
    return rcpp_result_gen;
END_RCPP
}
// estimateNetwork
Rcpp::List estimateNetwork(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::MatrixXd> M1, Eigen::Map<Eigen::MatrixXd> alpha, Eigen::Map<Eigen::MatrixXd> beta, Eigen::Map<Eigen::MatrixXd> gamma, Eigen::Map<Eigen::MatrixXd> tauAlpha, Eigen::Map<Eigen::MatrixXd> tauBeta, double rho, double lambda1a, double lambda1b, double lambda1g, double lambda2a, double lambda2b, Eigen::Map<Eigen::MatrixXd> XtX, Eigen::Map<Eigen::MatrixXd> XtXPlusRhoInv, Eigen::Map<Eigen::MatrixXd> XtM1, Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv, Eigen::Map<Eigen::MatrixXd> M1tY, Eigen::Map<Eigen::MatrixXd> XtY, Eigen::Map<Eigen::MatrixXd> L);
RcppExport SEXP _HDMAADMM_estimateNetwork(SEXP XSEXP, SEXP YSEXP, SEXP M1SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauAlphaSEXP, SEXP tauBetaSEXP, SEXP rhoSEXP, SEXP lambda1aSEXP, SEXP lambda1bSEXP, SEXP lambda1gSEXP, SEXP lambda2aSEXP, SEXP lambda2bSEXP, SEXP XtXSEXP, SEXP XtXPlusRhoInvSEXP, SEXP XtM1SEXP, SEXP M1tM1PlusRhoInvSEXP, SEXP M1tYSEXP, SEXP XtYSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauAlpha(tauAlphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauBeta(tauBetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1a(lambda1aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1b(lambda1bSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1g(lambda1gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2a(lambda2aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2b(lambda2bSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtXPlusRhoInv(XtXPlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtM1(XtM1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tM1PlusRhoInv(M1tM1PlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tY(M1tYSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtY(XtYSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateNetwork(X, Y, M1, alpha, beta, gamma, tauAlpha, tauBeta, rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b, XtX, XtXPlusRhoInv, XtM1, M1tM1PlusRhoInv, M1tY, XtY, L));
    return rcpp_result_gen;
END_RCPP
}
// estimatePathwayLasso
Rcpp::List estimatePathwayLasso(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::MatrixXd> M1, Eigen::Map<Eigen::MatrixXd> alpha, Eigen::Map<Eigen::MatrixXd> beta, Eigen::Map<Eigen::MatrixXd> gamma, Eigen::Map<Eigen::MatrixXd> tauAlpha, Eigen::Map<Eigen::MatrixXd> tauBeta, double rho, double lambda1a, double lambda1b, double lambda1g, double lambda2a, double lambda2b, Eigen::Map<Eigen::MatrixXd> XtX, Eigen::Map<Eigen::MatrixXd> XtXPlusRhoInv, Eigen::Map<Eigen::MatrixXd> XtM1, Eigen::Map<Eigen::MatrixXd> M1tM1PlusRhoInv, Eigen::Map<Eigen::MatrixXd> M1tY, Eigen::Map<Eigen::MatrixXd> XtY, double kappa, double nu);
RcppExport SEXP _HDMAADMM_estimatePathwayLasso(SEXP XSEXP, SEXP YSEXP, SEXP M1SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP tauAlphaSEXP, SEXP tauBetaSEXP, SEXP rhoSEXP, SEXP lambda1aSEXP, SEXP lambda1bSEXP, SEXP lambda1gSEXP, SEXP lambda2aSEXP, SEXP lambda2bSEXP, SEXP XtXSEXP, SEXP XtXPlusRhoInvSEXP, SEXP XtM1SEXP, SEXP M1tM1PlusRhoInvSEXP, SEXP M1tYSEXP, SEXP XtYSEXP, SEXP kappaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauAlpha(tauAlphaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type tauBeta(tauBetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1a(lambda1aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1b(lambda1bSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1g(lambda1gSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2a(lambda2aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda2b(lambda2bSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtXPlusRhoInv(XtXPlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtM1(XtM1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tM1PlusRhoInv(M1tM1PlusRhoInvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type M1tY(M1tYSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type XtY(XtYSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(estimatePathwayLasso(X, Y, M1, alpha, beta, gamma, tauAlpha, tauBeta, rho, lambda1a, lambda1b, lambda1g, lambda2a, lambda2b, XtX, XtXPlusRhoInv, XtM1, M1tM1PlusRhoInv, M1tY, XtY, kappa, nu));
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

static const R_CallMethodDef CallEntries[] = {
    {"_HDMAADMM_estimateElasticNet", (DL_FUNC) &_HDMAADMM_estimateElasticNet, 20},
    {"_HDMAADMM_estimateNetwork", (DL_FUNC) &_HDMAADMM_estimateNetwork, 21},
    {"_HDMAADMM_estimatePathwayLasso", (DL_FUNC) &_HDMAADMM_estimatePathwayLasso, 22},
    {"_HDMAADMM_fMatProd", (DL_FUNC) &_HDMAADMM_fMatProd, 3},
    {"_HDMAADMM_fMatTransProd", (DL_FUNC) &_HDMAADMM_fMatTransProd, 3},
    {"_HDMAADMM_fMatInv", (DL_FUNC) &_HDMAADMM_fMatInv, 2},
    {"_HDMAADMM_fMatChol", (DL_FUNC) &_HDMAADMM_fMatChol, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_HDMAADMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
