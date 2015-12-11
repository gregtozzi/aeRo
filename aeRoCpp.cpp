#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List naca4 (int m100, int p10, int t100, int panels, bool hcos) {
  
  // Scale the input
  double m = m100 / 100.0;
  double p = p10 / 10.0;
  double t = t100 / 100.0;
  
  // Establish the x-location of the verticies
  arma::vec x;
  int numVert = panels + 1;
  if (hcos == TRUE) {
    arma::vec theta = linspace(0, datum::pi / 2, numVert);
    x = - cos(theta) + 1;
  } else {
    x = linspace(0, 1, numVert);
  }
  
  // Build the camber line
  arma::vec zc(numVert, fill::zeros);
  arma::vec dzcdx(numVert, fill::zeros);
  for (int i = 0; i < numVert; i++) {
    if (x(i) <= p) {
      zc(i) = m * x(i) / pow(p, 2) * (2 * p - x(i));
      dzcdx(i) = 2 * m / pow(p, 2) * (p - x(i));
    } else {
      zc(i) = m * (1 - x(i)) / pow(1 - p, 2) * (1 + x(i) - 2 * p);
      dzcdx(i) = 2 * m / pow(1 - p, 2) * (p - x(i));
    }
  }
  
  // Apply the thickness
  arma::vec theta = atan(dzcdx);
  arma::vec zt = 5 * t * (0.2969 * sqrt(x)
                            - 0.126 * x
                            - 0.3516 * pow(x, 2)
                            + 0.2843 * pow(x, 3)
                            - 0.1036 * pow(x, 4));
  arma::vec xu = x - zt % sin(theta);
  arma::vec xl = x + zt % sin(theta);
  arma::vec zu = zc + zt % cos(theta);
  arma::vec zl = zc - zt % cos(theta);
  zu(panels) = 0;
  zl(panels) = 0;
                            
  // Build the output
  List ret;
  ret["x"] = x;
  ret["zt"] = zt;
  ret["zc"] = zc;
  ret["xu"] = xu;
  ret["xl"] = xl;
  ret["zu"] = zu;
  ret["zl"] = zl;
  return(ret);
}