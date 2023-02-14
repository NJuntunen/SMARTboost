#include <Rcpp.h>

using namespace Rcpp;

// define lnpmu C++ function
double lnpmu(double mu, double varmu, double dofmu) {
  double s = sqrt(varmu);
  double lnp = R::dnorm(mu, 0.0, s, true);
  return lnp;
}

// define lnptau C++ function
double lnptau(double tau, double meanlntau, double varlntau, double doflntau, double depth) {
  double s = sqrt(varlntau/depth);
  double lnp = R::dt(log(tau), doflntau, true) - log(tau) + R::dnorm(log(tau), meanlntau, s, true);
  return lnp;
}

// define the main function
// [[Rcpp::export]]
List fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, List param, List infeaturesfit, bool dichotomous, double mu, double tau, bool dichotomous_i) {
  Eigen::VectorXd diagGGh = G.colwise().squaredNorm();
  Eigen::MatrixXd GGh = G.transpose() * G;
  int n = G.rows();
  int p = G.cols();
  double var_r = var_epsilon/(1-param["R2p"]);
  double Pb = diagGGh.sum()/(n*var_r*param["R2p"]);
  Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(p,p);
  Eigen::MatrixXd GGh_var_r_Pb_I_p = GGh + var_r*param["loglikdivide"]*Pb*I_p;
  Eigen::VectorXd beta(p);
  List ret;

  try {
    beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
  } catch (...) {
    while (true) {
      Pb = Pb * 2.01;
      GGh_var_r_Pb_I_p = GGh + var_r*param["loglikdivide"]*Pb*I_p;
      try {
        beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
        break;
      } catch (...) {
        // Do nothing, loop again
      }
    }
  }

  Eigen::VectorXd Gbeta = G * beta;
  double loglik = -0.5*(r-Gbeta).squaredNorm()/var_r/param["loglikdivide"];
  double logpdfbeta = -0.5*(p*log(2*M_PI) - p*log(Pb) + Pb*beta.squaredNorm());

  double logpdfmu = 0;
  double logpdftau = 0;

  if (!dichotomous_i) {
    if (param["sharptree"]) {
      logpdfmu = lnpmu(mu, param["varmu"], param["dofmu"]);
    } else {
      logpdfmu = lnpmu(mu, param["varmu"], param["dofmu"]);
      logpdftau = lnptau(tau, param["meanlntau"], param["varlntau"], param["doflntau"], param["depth"]);
    }
  }

  double loss = -(loglik + logpdfbeta + logpdftau + logpdfmu);

  ret["loss"] = loss;
  ret["Gbeta"] = Gbeta;
  ret["beta"] = beta;

  return ret;
}
