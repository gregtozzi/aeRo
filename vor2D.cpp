#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::rowvec vor2D (double G, double x, double z, double xj, double zj) {
  
  double xmxj = x - xj;
  double zmzj = z - zj;
  double rj2 = pow(xmxj, 2) + pow(zmzj, 2);
  double premult = G / 2. / datum::pi / rj2;
  double u = premult * zmzj;
  double w = premult * (-1) * xmxj;
  
  arma::rowvec ret;
  ret << u << w;
  return(ret);
}

// [[Rcpp::export]]
List lumpVortex (arma::vec xc, arma::vec zc,
                 arma::vec x, arma::vec z,
                 arma::vec nx, arma::vec nz,
                 float Q, float alpha) {
  
  // Initialize output and intermediate matricies
  int matDim = x.n_elem;
  arma::mat u(matDim, matDim);
  arma::mat w(matDim, matDim);
  arma::mat a(matDim, matDim);
  
  // Initialize intermediate vectors
  arma::rowvec vorOut(2);
  arma::colvec N(2);
  arma::colvec RHS(matDim);
  arma::vec freeStream(2);
  arma::colvec G(matDim);
  
  // Decompose the free stream velocity
  freeStream << - Q * cos(alpha) << - Q * sin(alpha);
  
  // Initialize iterators
  int i;
  int j;
  
  // Compute the influence coefficients and the RHS
  for(i = 0; i < matDim; i++) {
    N << nx(i) << nz(i);
    RHS(i) = dot(freeStream, N);
    for(j = 0; j < matDim; j++) {
      vorOut = vor2D(1, xc(i), zc(i), x(j), z(j));
      a(i, j) = dot(vorOut, N);
    }
  }
  
  // Prepare the output
  List ret;
  ret["a"] = a;
  ret["RHS"] = RHS;
  ret["G"] = solve(a, RHS); // Solve for the vortex strengths
  return(ret);
}

// [[Rcpp::export]]
arma::mat vor2DV (double G, arma::vec x, arma::vec z, arma::vec xj, arma::vec zj) {
  
  arma::vec xmxj = x - xj;
  arma::vec zmzj = z - zj;
  arma::vec rj2 = pow(xmxj, 2) + pow(zmzj, 2);
  arma::vec premult = G / 2. / datum::pi / rj2;
  arma::vec u = premult % zmzj;
  arma::vec w = premult * (-1) % xmxj;
  
  arma::mat ret(x.n_elem, 2);
  ret.col(0) = u;
  ret.col(1) = w;
  return(ret);
}

// [[Rcpp::export]]
List vor2DM (double G, arma::mat x, arma::mat z, double xj, double zj) {
  
  arma::mat xmxj = x - xj;
  arma::mat zmzj = z - zj;
  arma::mat rj2 = pow(xmxj, 2) + pow(zmzj, 2);
  arma::mat premult = G / 2. / datum::pi / rj2;
  arma::mat u = premult % zmzj;
  arma::mat w = premult * (-1) % xmxj;
  
  List ret;
  ret["u"] = u;
  ret["w"] = w;
  return(ret);
}