// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// Fast crossproduct of single matrix transpose(X)*X
// [[Rcpp::export]]
Eigen::MatrixXf fcprd(const Eigen::MatrixXf X){
    const int n = X.cols();
    return Eigen::MatrixXf(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}
// Fast crossproduct of two matrices
// [[Rcpp::export]]
Eigen::MatrixXf cpp_prod(const Eigen::MatrixXf X, const Eigen::MatrixXf Y){
  return Eigen::MatrixXf(X*Y);
}
// Cholesky decomposition of a PD matrix (outputs upper triangular matrix)
// [[Rcpp::export]]
Eigen::MatrixXf cpp_cholesky(const Eigen::MatrixXf M){
  return Eigen::MatrixXf(M.llt().matrixU());
}