#include <Rcpp.h>
#include <Eigen/Core>
#include <RcppEigen.h>
#include <iostream>

using namespace Rcpp;

// sigmoidf function
// [[Rcpp::export]]
NumericVector sigmoidf_cpp(NumericVector x, double mu, double tau, std::string sigmoid, bool dichotomous) {

  NumericVector g ;
  if(dichotomous){
    for (int i = 0; i < x.size(); i++) {
      g(i) = (x(i) > 0) ? 1 : 0;
    }
    return g;
  } else {
    if(sigmoid == "sigmoidsqrt"){
      g = 0.5 + 0.5 * ( 0.5 * tau * (x - mu) / sqrt((1 + pow(0.5 * tau * (x - mu),2))));
      // Eigen::VectorXd g = 0.5 + 0.5 * (0.5 * tau * (x.array() - mu) / sqrt(1 + (0.5 * tau * (x.array() - mu))*(0.5 * tau * (x.array() - mu))));
    } else if(sigmoid == "sigmoidlogistic"){
      g = 1 - 1/(1 + exp(tau * (x - mu)));
    }
  }

  return g;
}

// updateG_allocated function
// [[Rcpp::export]]
Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, NumericVector g, Eigen::MatrixXd G) {
  int n = G0.rows(), p = G0.cols();
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < n; j++) {
      G(j, i) = G0(j, i) * g[j];
      G(j, i + p) = G0(j, i) * (1 - g[j]);
    }
  }
  return G;
}

// define lnpmu C++ function
double lnpmu(double mu, double varmu, double dofmu) {
  double s = sqrt(varmu);
  double lnp = R::dnorm(mu, 0.0, s, true);
  return lnp;
}

// logpdft function
// x: input values
// m: mean
// s: standard deviation
// v: degrees of freedom
NumericVector logpdft(NumericVector x, double m, double s, double v) {
  double constant = -0.5723649429247001 + lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v);
  NumericVector z = (x - m) / s;
  NumericVector logpdfz = constant - 0.5 * (1 + v) * log(1 + pow(z, 2) / v);
  return logpdfz - log(s);
}

// lnptau function
// tau: input values
// meanlntau: mean of log(tau)
// varlntau: variance of log(tau)
// doflntau: degrees of freedom
// depth: tree depth
double lnptau(double tau, double meanlntau, double varlntau, double doflntau, double depth) {
  double s = sqrt(varlntau / depth);
  NumericVector logpdf_tau = logpdft(log(tau), meanlntau, s, doflntau);
  double lnp_tau = sum(logpdf_tau);
  return lnp_tau;
}


// define the main function
// [[Rcpp::export]]
List fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, List param, List infeaturesfit,
                 List dichotomous, double mu, double tau, bool dichotomous_i) {

  Eigen::VectorXd diagGGh = G.colwise().squaredNorm();
  Eigen::MatrixXd GGh = G.transpose() * G;
  int n = G.rows();
  int p = G.cols();
  double var_r = var_epsilon/(1-as<double>(param["R2p"]));
  double Pb = diagGGh.sum()/(n*var_r*as<double>(param["R2p"]));
  Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(p,p);
  Eigen::MatrixXd GGh_var_r_Pb_I_p = GGh + var_r*as<double>(param["loglikdivide"])*Pb*I_p;
  Eigen::VectorXd beta(p);
  List ret;

  try {
    beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
  } catch (...) {
    while (true) {
      Pb = Pb * 2.01;
      GGh_var_r_Pb_I_p = GGh + var_r*as<double>(param["loglikdivide"])*Pb*I_p;
      try {
        beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
        break;
      } catch (...) {
        // Do nothing, loop again
      }
    }
  }

  Eigen::VectorXd Gbeta = G * beta;
  double loglik = -0.5*(r-Gbeta).squaredNorm()/var_r/as<double>(param["loglikdivide"]);
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



// Gfitbeta function
// [[Rcpp::export]]
double Gfitbeta_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, NumericVector xi, List param,
                    double var_epsilon, List infeaturesfit, List dichotomous, List muLogTau,
                    bool dichotomous_i, Eigen::MatrixXd G) {
  NumericVector gL;
  double mu = as<double>(muLogTau[0]);
  double tau = exp(as<double>(muLogTau[1]));
  tau = std::max(tau, 0.2);

  gL = sigmoidf_cpp(xi, mu, tau, param["sigmoid"], dichotomous_i);

  G = updateG_allocated_cpp(G0, gL, G);

  List result = fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit,
                            dichotomous, mu, tau, dichotomous_i);

  return result[0];
}















