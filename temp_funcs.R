// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <Rcpp.h>
#include <Eigen/Core>
#include <RcppEigen.h>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <RcppThread.h>


using namespace Rcpp;
using namespace RcppThread;

// sigmoidf function
// [[Rcpp::export]]
NumericVector sigmoidf_cpp(NumericVector x, double mu, double tau, std::string sigmoid, bool dichotomous) {

  NumericVector g(x.size());
  if(dichotomous){
    for (int i = 0; i < x.size(); i++) {
      g(i) = (x(i) > 0) ? 1 : 0;
    }
    return g;
  } else {
    if(sigmoid == "sigmoidsqrt"){
      g = 0.5 + 0.5 * ( 0.5 * tau * (x - mu) / sqrt((1 + pow(0.5 * tau * (x - mu),2))));
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
NumericVector logpdft(Eigen::VectorXd x, double m, double s, double v) {
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
double lnptau(double tau, double meanlntau, double varlntau, double doflntau, int depth) {
  double s = sqrt(varlntau / depth);
  NumericVector logpdf_tau = logpdft(log(tau), meanlntau, s, doflntau);
  double lnp_tau = sum(logpdf_tau);
  return lnp_tau;
}


// define the main function
// [[Rcpp::export]]
List fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, List param, LogicalVector infeaturesfit,
                 LogicalVector dichotomous, double mu, double tau, bool dichotomous_i,
                 double R2p, double loglikdivide, bool sharptree, double varmu, double varlntau, double meanlntau, double dofmu, double doflntau,  int depth) {

  Eigen::VectorXd diagGGh = G.colwise().squaredNorm();
  Eigen::MatrixXd GGh = G.transpose() * G;
  int n = G.rows();
  int p = G.cols();
  double var_r = var_epsilon/(1-R2p);
  double Pb = diagGGh.sum()/(n*var_r*R2p);
  Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(p,p);
  Eigen::MatrixXd GGh_var_r_Pb_I_p = GGh + var_r*loglikdivide*Pb*I_p;
  Eigen::VectorXd beta(p);
  List ret;

  try {
    beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
  } catch (...) {
    while (true) {
      Pb = Pb * 2.01;
      GGh_var_r_Pb_I_p = GGh + var_r*loglikdivide*Pb*I_p;
      try {
        beta = GGh_var_r_Pb_I_p.llt().solve(G.transpose() * r);
        break;
      } catch (...) {
        // Do nothing, loop again
      }
    }
  }

  Eigen::VectorXd Gbeta = G * beta;
  double loglik = -0.5*(r-Gbeta).squaredNorm()/var_r/loglikdivide;
  double logpdfbeta = -0.5*(p*log(2*M_PI) - p*log(Pb) + Pb*beta.squaredNorm());

  double logpdfmu = 0;
  double logpdftau = 0;

  if (!dichotomous_i) {
    if (sharptree) {
      logpdfmu = lnpmu(mu, varmu, dofmu);
    } else {
      logpdfmu = lnpmu(mu, varmu, dofmu);
      logpdftau = lnptau(tau, meanlntau, varlntau, doflntau, depth);
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
                    double var_epsilon, LogicalVector infeaturesfit, LogicalVector dichotomous, List muLogTau,
                    bool dichotomous_i, Eigen::MatrixXd G,
                    double R2p, double loglikdivide, bool sharptree, double varmu, double dofmu, double varlntau, double meanlntau, double doflntau,  int depth) {
  NumericVector gL;
  double mu = as<double>(muLogTau[0]);
  double tau = exp(as<double>(muLogTau[1]));
  tau = std::max(tau, 0.2);

  gL = sigmoidf_cpp(xi, mu, tau, param["sigmoid"], dichotomous_i);

  G = updateG_allocated_cpp(G0, gL, G);

  List result = fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit,
                            dichotomous, mu, tau, dichotomous_i,
                            R2p, loglikdivide, sharptree, varmu,varlntau, meanlntau, dofmu, doflntau, depth);

  return result[0];
}


// [[Rcpp::export]]
DataFrame add_depth_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, NumericVector xi, LogicalVector infeaturesfit,
                        LogicalVector dichotomous, bool dichotomous_i, NumericVector mugridi, NumericVector taugrid, List param, double var_epsilon,
                        double R2p, double loglikdivide, bool sharptree, double varmu,double varlntau, double meanlntau, double dofmu, double doflntau,  int depth) {


  NumericMatrix lossmatrix(taugrid.length(), mugridi.length());
  std::fill(lossmatrix.begin(), lossmatrix.end(), R_PosInf);
  int n = G0.rows();
  int p = G0.cols();
  Eigen::MatrixXd G(n, 2*p);

  double loss;
  double tau;
  double mu;

  if(dichotomous_i) {
    // no optimization needed
    loss = Gfitbeta_cpp(r, h, G0, xi, param, var_epsilon, infeaturesfit, dichotomous, List::create(0,0), dichotomous_i, G,
                        R2p, loglikdivide, sharptree, varmu,varlntau, meanlntau, dofmu, doflntau, depth);
    tau = 999.9;
    mu = 0;
  } else {
    for(int indexmu = 0; indexmu < mugridi.size(); indexmu++) {
      for(int indextau = 0; indextau < taugrid.size(); indextau++) {
        lossmatrix(indextau, indexmu) = Gfitbeta_cpp(r, h, G0, xi, param, var_epsilon, infeaturesfit, dichotomous, List::create(mugridi[indexmu],log(taugrid[indextau])), dichotomous_i, G,
                                                     R2p, loglikdivide, sharptree, varmu, varlntau, meanlntau, dofmu, doflntau, depth);
        if (indextau > 0 && lossmatrix(indextau, indexmu) > lossmatrix(indextau-1, indexmu)) {
          break;
        }
      }
    }

    int min_row = 0;
    int min_col = 0;
    double min_val = lossmatrix(0, 0);

    for (int i = 0; i < lossmatrix.nrow(); i++) {
      for (int j = 0; j < lossmatrix.ncol(); j++) {
        double val = lossmatrix(i, j);
        if (val < min_val) {
          min_val = val;
          min_row = i;
          min_col = j;
        }
      }
    }

    loss = lossmatrix(min_row, min_col);
    tau = taugrid(min_row);
    mu = mugridi(min_col);

    // Optionally, further optimize over mu
    if(param["optimizevs"]) {

      //define nlopt optimization here

    }
  }

  return DataFrame::create(_["loss"] = loss, _["tau"] = tau, _["mu"] = mu);
}

LogicalVector updateinfeatures_cpp(LogicalVector infeatures, NumericVector ifit) {
  LogicalVector x = clone(infeatures);
  for (auto &i : ifit) {
    x[i] = TRUE;
  }
  return x;
}

// Worker function to process each column
// 'outputarray' is locked before accessing/modifying it
void loop_column(int i, int ps, NumericMatrix& outputarray,
                 const Eigen::VectorXd& r, const Eigen::VectorXd& h, const Eigen::MatrixXd& G0,
                 const NumericMatrix& x, const LogicalVector& infeatures, const LogicalVector& dichotomous,
                 const NumericMatrix& mugrid, const NumericVector& taugrid, const List& param, double var_epsilon,
                 double R2p, double loglikdivide, bool sharptree, double varmu, double varlntau, double meanlntau, double dofmu, double doflntau,  int depth) {

  NumericVector xi = x(_, ps);
  LogicalVector infeaturesfit = updateinfeatures_cpp(infeatures, ps);
  bool dichotomous_i = dichotomous[ps];
  NumericVector mugridi = mugrid(_, ps);

  DataFrame output = add_depth_cpp(r, h, G0, xi, infeaturesfit,
                                   dichotomous, dichotomous_i, mugridi, taugrid, param, var_epsilon,
                                   R2p, loglikdivide, sharptree, varmu, dofmu, varlntau, meanlntau, doflntau, depth);

  outputarray(i, 0) = output["loss"];
  outputarray(i, 1) = output["tau"];
  outputarray(i, 2) = output["mu"];
}

// [[Rcpp::export]]
NumericMatrix loopfeatures_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0,
                               NumericMatrix x, NumericVector ifit, LogicalVector infeatures,
                               NumericMatrix mugrid, LogicalVector dichotomous, NumericVector taugrid,
                               List param, double var_epsilon) {
  int p = x.ncol();
  NumericMatrix outputarray(p, 3); // [loss, τ, μ]

  double R2p = param["R2p"];
  double loglikdivide = param["loglikdivide"];
  bool sharptree = param["sharptree"];
  double varmu = param["varmu"];
  double meanlntau = param["meanlntau"];
  double varlntau = param["varlntau"];
  double dofmu = param["dofmu"];
  double doflntau = param["doflntau"];
  int depth = param["depth"];

  IntegerVector ps(p);
  for (int i = 0; i < p; i++) {
    ps[i] = i+1;
  }

  if (as<double>(param["subsampleshare_columns"]) < 1) {
    int psmall = round(p * as<double>(param["subsampleshare_columns"]));
    IntegerVector p_sampled = Rcpp::sample(ps, psmall);
    ps = p_sampled;
  }

  // loop through the columns in parallel
  RcppThread::parallelFor(0, p, [&](int i) {
    loop_column(i,ps[i]-1, outputarray, r, h, G0, x, infeatures, dichotomous, mugrid, taugrid, param, var_epsilon,
                R2p, loglikdivide, sharptree, varmu, meanlntau, varlntau, dofmu, doflntau, depth);
  });

  return outputarray;
}


// Worker function to process each column
// 'outputarray' is locked before accessing/modifying it
void loop_column(int i, int ps, NumericMatrix& outputarray,
                 const Eigen::VectorXd& r, const Eigen::VectorXd& h, const Eigen::MatrixXd& G0,
                 const Eigen::VectorXd& x, const LogicalVector& infeatures, const LogicalVector& dichotomous,
                 const Eigen::MatrixXd& mugrid, const NumericVector& taugrid, const List& param, double var_epsilon) {

  try {
    Eigen::VectorXd xi = x.col(ps);
    LogicalVector infeaturesfit = updateinfeatures_cpp(infeatures, ps);
    bool dichotomous_i = dichotomous[ps];
    Eigen::VectorXd mugridi = mugrid.col(ps);

    // DataFrame output = add_depth_cpp(r, h, G0, xi, infeaturesfit,
                                        //                                  dichotomous, dichotomous_i, mugridi, taugrid, param, var_epsilon);
    double loss;
    double tau;
    double mu;

    loss= 100;
    tau = 999.9;
    mu = 0;


    DataFrame output =  DataFrame::create(_["loss"] = loss, _["tau"] = tau, _["mu"] = mu);

    outputarray(i, 0) = output["loss"];
    outputarray(i, 1) = output["tau"];
    outputarray(i, 2) = output["mu"];

  } catch (std::exception& e) {
    Rcpp::stop("Error in loop_column: %s", e.what());
  }
}

