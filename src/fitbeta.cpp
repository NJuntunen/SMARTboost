
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <RcppEigen.h>
#include <RcppThread.h>
#include <cmath>
#include "structures.h"


using namespace Rcpp;
// using namespace RcppParallel;

// sigmoidf function
Eigen::VectorXd sigmoidf_cpp(Eigen::VectorXd x, double mu, double tau, std::string sigmoid, bool dichotomous) {

  Eigen::VectorXd g(x.size());
  if(dichotomous){
    for (int i = 0; i < x.size(); i++) {
      g(i) = (x(i) > 0) ? 1 : 0;
    }
    return g;
  } else {
    if(sigmoid == "sigmoidsqrt"){
      g = 0.5 + 0.5 * ( 0.5 * tau * (x.array() - mu) / sqrt((1 + pow(0.5 * tau * (x.array() - mu),2))));
    } else if(sigmoid == "sigmoidlogistic"){
      g = 1 - 1/(1 + exp(tau * (x.array() - mu)));
    }
  }

  return g;
}

// updateG_allocated function

Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, Eigen::VectorXd g, Eigen::MatrixXd G) {

  int p = G0.cols();
  for (int i = 0; i < p; i++) {
    G.col(i) = G0.col(i).array() * g.array();
    G.col(i + (p-1)) = G0.col(i).array() * (1 - g.array());
  }
  // int n = G0.rows(), p = G0.cols();
  // for (int i = 0; i < p; i++) {
  //   for (int j = 0; j < n; j++) {
  //     G(j, i) = G0(j, i) * g[j];
  //     G(j, i + p) = G0(j, i) * (1 - g[j]);
  //   }
  // }
  return G;
}

// define lnpmu C++ function
double lnpmu(double mu, double varmu, double dofmu) {
  double s = sqrt(varmu);
  double lnp = exp(-0.5*pow((mu-0)/s,2))/(sqrt(2*M_PI)*s);
  return lnp;
}

// logpdft function
// x: input values
// m: mean
// s: standard deviation
// v: degrees of freedom
double logpdft(double x, double m, double s, double v) {
  double constant = -0.5723649429247001 + lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v);
  double z = (x - m) / s;
  double logpdfz = constant - 0.5 * (1 + v) * log(1 + pow(z, 2) / v);
  return logpdfz - log(s);
}
//
// lnptau function
// tau: input values
// meanlntau: mean of log(tau)
// varlntau: variance of log(tau)
// doflntau: degrees of freedom
// depth: tree depth
double lnptau(double tau, double meanlntau, double varlntau, double doflntau, int depth) {
  double s = sqrt(varlntau / depth);
  double logpdf_tau = logpdft(log(tau), meanlntau, s, doflntau);
  // double lnp_tau = logpdf_tau;
  return logpdf_tau;
}


// define the main function
FitBetaStruct fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, SMARTParamStruct SMARTparams, double mu, double tau, bool dichotomous_i) {

  Eigen::VectorXd diagGGh = G.colwise().squaredNorm();
  Eigen::MatrixXd GGh = G.transpose() * G;
  int n = G.rows();
  int p = G.cols();
  double var_r = var_epsilon/(1-SMARTparams.R2p);
  double Pb = diagGGh.sum()/(n*var_r*SMARTparams.R2p);
  Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(p,p);
  Eigen::MatrixXd GGh_var_r_Pb_I_p = GGh + var_r*SMARTparams.loglikdivide*Pb*I_p;
  Eigen::VectorXd beta(p);

  // Use BiCGSTAB to solve the system
  Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
  solver.compute(GGh_var_r_Pb_I_p);
  beta = solver.solve(G.transpose() * r);


  // Check if the solution is valid, and if not, use a larger Pb
  // int max_iter = 5; // set a maximum number of iterations
  // int iter = 0;
  while (solver.info() != Eigen::Success) {
    // iter++;
    Pb = Pb * 2.01;
    GGh_var_r_Pb_I_p = GGh + var_r * SMARTparams.loglikdivide * Pb * I_p;
    solver.compute(GGh_var_r_Pb_I_p);
    beta = solver.solve(G.transpose() * r);
  }

  Eigen::VectorXd Gbeta = G * beta;
  double loglik = -0.5*(r-Gbeta).squaredNorm()/var_r/SMARTparams.loglikdivide;
  double logpdfbeta = -0.5*(p*log(2*M_PI) - p*log(Pb) + Pb*beta.squaredNorm());

  double logpdfmu = 0;
  double logpdftau = 0;

  if (!dichotomous_i) {
    if (SMARTparams.sharptree) {
      logpdfmu = lnpmu(mu, SMARTparams.varmu, SMARTparams.dofmu);
    } else {
      logpdfmu = lnpmu(mu, SMARTparams.varmu, SMARTparams.dofmu);
      logpdftau = lnptau(tau, SMARTparams.meanlntau, SMARTparams.varlntau, SMARTparams.doflntau, SMARTparams.depth);
    }
  }

  double loss = -(loglik + logpdfbeta + logpdftau + logpdfmu);

  FitBetaStruct FitBeta;
  FitBeta.loss = loss;
  FitBeta.Gbeta = Gbeta;
  FitBeta.beta = beta;

  return FitBeta;
}



// Gfitbeta function
double Gfitbeta_cpp(const Eigen::VectorXd r, const Eigen::VectorXd h, const Eigen::MatrixXd G0, const Eigen::VectorXd xi, SMARTParamStruct SMARTparams,
                    const double var_epsilon, std::vector<double> muLogTau,
                    const bool dichotomous_i, Eigen::MatrixXd G) {
  Eigen::VectorXd gL(xi.size());
  double tau = exp(muLogTau[1]);
  tau = std::max(tau, 0.2);


  gL = sigmoidf_cpp(xi, muLogTau[0], tau, SMARTparams.sigmoid, dichotomous_i);

  G = updateG_allocated_cpp(G0, gL, G);

  FitBetaStruct result = fitbeta_cpp(r, G, var_epsilon, SMARTparams, muLogTau[0], tau, dichotomous_i);

  return result.loss;
}



// Function that calculates the loss, tau, and mu for one column of x
std::vector<double> add_depth_cpp(const Eigen::VectorXd x, const Eigen::VectorXd r, const Eigen::VectorXd h,
                                  const Eigen::MatrixXd G0, const bool dichotomous_i, const Eigen::VectorXd mugridi,
                                  const Eigen::VectorXd taugrid, SMARTParamStruct SMARTparams, const double var_epsilon) {




    // create a 2D vector of size nrow x ncol with default value R_PosInf

  const int n = G0.rows();
  const int p = G0.cols();

  std::vector<double> result(3);
  Eigen::MatrixXd lossmatrix(taugrid.size(), mugridi.size());
  Eigen::MatrixXd G(n, 2*p);
  std::vector<double> muLogTau(2);

  lossmatrix.fill(std::numeric_limits<double>::infinity());

  double loss;
  double tau;
  double mu;

  if(dichotomous_i) {

    muLogTau[0] = 0;
    muLogTau[1] = 0;
    // no optimization needed
    loss = Gfitbeta_cpp(r, h, G0, x, SMARTparams, var_epsilon, muLogTau, dichotomous_i, G);
    tau = 999.9;
    mu = 0;
  } else {
    for(int indexmu = 0; indexmu < mugridi.size(); indexmu++) {
      for(int indextau = 0; indextau < taugrid.size(); indextau++) {

        muLogTau[0] = mugridi[indexmu];
        muLogTau[1] = log(taugrid[indextau]);

        lossmatrix(indextau, indexmu) = Gfitbeta_cpp(r, h, G0, x, SMARTparams, var_epsilon, muLogTau, dichotomous_i, G);
        if (indextau > 0 && lossmatrix(indextau, indexmu) > lossmatrix(indextau-1, indexmu)) {
          break;
        }
      }
    }

    int min_row = 0;
    int min_col = 0;
    double min_val = lossmatrix(0, 0);

    for (int i = 0; i < lossmatrix.rows(); i++) {
      for (int j = 0; j < lossmatrix.cols(); j++) {
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
    if(SMARTparams.optimizevs) {

      //define nlopt optimization here

    }
  }

  result[0] = loss;
  result[1] = tau;
  result[2] = mu;

  return result;
}

Eigen::MatrixXd loopfeatures_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0,
                        Eigen::MatrixXd x, Eigen::MatrixXd mugrid, Eigen::VectorXd dichotomous, Eigen::VectorXd taugrid,
                        SMARTParamStruct SMARTparams, double var_epsilon) {

  int p = x.cols();
//   // std::cout << "The value of x is: " << p << std::endl;
  IntegerVector ps(p);

  for (int i = 0; i < p; i++) {
    ps[i] = i+1;
  }

  if (SMARTparams.subsampleshare_columns < 1) {
    int psmall = round(p * SMARTparams.subsampleshare_columns);
    IntegerVector p_sampled = Rcpp::sample(ps, psmall);
    ps = p_sampled;
  }

  // loop through the columns in parallel

  int count = 0;
  Eigen::MatrixXd output(ps.size(), 3);
  std::mutex mtx;

  RcppThread::parallelForEach(ps, [&](int i) {
    int index = i-1;
    Eigen::VectorXd xi = x.col(index);
    bool dichotomous_i = dichotomous[index];
    Eigen::VectorXd mugridi = mugrid.col(index);

    std::vector<double> result = add_depth_cpp(xi, r, h, G0,
                                               dichotomous_i, mugridi, taugrid, SMARTparams, var_epsilon);

    mtx.lock();
    int current_count = count;
    output(current_count, 0) = result[0];
    output(current_count, 1) = result[1];
    output(current_count, 2) = result[2];
    count = count+1;
    mtx.unlock();
  });

  RcppThread::wait();

  return output;

}






