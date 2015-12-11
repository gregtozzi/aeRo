#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec vor2D (float G, double x, double z, double xj, double zj) {
  
  double xmxj = x - xj;
  double zmzj = z - zj;
  double rj2 = pow(xmxj, 2) + pow(zmzj, 2);
  double premult = G / 2. / datum::pi / rj2;
  double u = premult * zmzj;
  double w = premult * (-1) * xmxj;
  
  arma::vec ret;
  ret << u << w;
  return(ret);
}

// [[Rcpp::export]]
arma::vec wrapper (float G, double x, double z, double xj, double zj) {
  arma::vec y = vor2D(G, x, z, xj, zj);
  return(y);
}