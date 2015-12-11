#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List joukowsky (int cenChi, int cenEta, int radius, int points) {
  
  // Build the circle in the zeta plane
  float lastPoint = 2 * datum::pi * (1 - 1 / points); 
  arma::vec theta;
  theta = linspace(0, lastPoint, points);
  arma::vec chi = sin(theta) * radius + cenChi;
  arma::vec eta = cos(theta) * radius + cenEta;
  
  // Map the circle into the z plane
  arma::vec x = chi % (pow(chi, 2) + pow(eta, 2) + 1) / (pow(chi, 2) + pow(eta, 2));
  arma::vec z = eta % (pow(chi, 2) + pow(eta, 2) - 1) / (pow(chi, 2) + pow(eta, 2));
  
  // Return the output
  List ret;
  ret["x"] = x;
  ret["z"] = z;
  return(ret);
}