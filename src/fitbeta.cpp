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


// sigmoidf function
// [[Rcpp::export]]
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
// [[Rcpp::export]]
Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, Eigen::VectorXd g, Eigen::MatrixXd G) {
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
double lnptau(double tau, double meanlntau, double varlntau, double doflntau, int depth) {
  double s = sqrt(varlntau / depth);
  NumericVector logpdf_tau = logpdft(log(tau), meanlntau, s, doflntau);
  double lnp_tau = sum(logpdf_tau);
  return lnp_tau;
}


// define the main function
// [[Rcpp::export]]
List fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, List param, LogicalVector infeaturesfit,
                 double mu, double tau, bool dichotomous_i) {

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

  // Use BiCGSTAB to solve the system
  Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
  solver.compute(GGh_var_r_Pb_I_p);
  beta = solver.solve(G.transpose() * r);

  // Check if the solution is valid, and if not, use a larger Pb
  if (solver.info() != Eigen::Success) {
    while (true) {
      Pb = Pb * 2.01;
      GGh_var_r_Pb_I_p = GGh + var_r*as<double>(param["loglikdivide"])*Pb*I_p;
      solver.compute(GGh_var_r_Pb_I_p);
      beta = solver.solve(G.transpose() * r);
      if (solver.info() == Eigen::Success) {
        break;
      }
      // Do nothing, loop again
    }
  }

  Eigen::VectorXd Gbeta = G * beta;
  double loglik = -0.5*(r-Gbeta).squaredNorm()/var_r/as<double>(param["loglikdivide"]);
  double logpdfbeta = -0.5*(p*log(2*M_PI) - p*log(Pb) + Pb*beta.squaredNorm());

  double logpdfmu = 0;
  double logpdftau = 0;

  if (!dichotomous_i) {
    if (as<double>(param["sharptree"])) {
      logpdfmu = lnpmu(mu, as<double>(param["varmu"]), as<double>(param["dofmu"]));
    } else {
      logpdfmu = lnpmu(mu, as<double>(param["varmu"]), as<double>(param["dofmu"]));
      logpdftau = lnptau(tau, as<double>(param["meanlntau"]), as<double>(param["varlntau"]), as<double>(param["doflntau"]), as<double>(param["depth"]));
    }
  }

  double loss = -(loglik + logpdfbeta + logpdftau + logpdfmu);

  ret["loss"] = loss;
  ret["Gbeta"] = Gbeta;
  ret["beta"] = beta;

  return ret;
}



// Gfitbeta function
double Gfitbeta_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi, List param,
                    double var_epsilon, LogicalVector infeaturesfit, List muLogTau,
                    bool dichotomous_i, Eigen::MatrixXd G) {
  Eigen::VectorXd gL(xi.size());
  double tau = exp(as<double>(muLogTau[1]));
  tau = std::max(tau, 0.2);


  gL = sigmoidf_cpp(xi, as<double>(muLogTau[0]), tau, param["sigmoid"], dichotomous_i);

  G = updateG_allocated_cpp(G0, gL, G);

  List result = fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit, as<double>(muLogTau[0]), tau, dichotomous_i);


  return 10;
}

LogicalVector updateinfeatures_cpp(LogicalVector infeatures, NumericVector ifit) {
  LogicalVector x = clone(infeatures);
  for (auto &i : ifit) {
    x[i] = TRUE;
  }
  return x;
}

// Function that calculates the loss, tau, and mu for one column of x
Eigen::VectorXd add_depth_cpp(Eigen::VectorXd& x, Eigen::VectorXd& r, Eigen::VectorXd& h,
                              Eigen::MatrixXd& G0, LogicalVector& infeaturesfit,
                              bool& dichotomous_i, Eigen::VectorXd& mugridi,
                              NumericVector& taugrid, List& param, double& var_epsilon) {

  Eigen::VectorXd result(3);

  NumericMatrix lossmatrix(taugrid.length(), mugridi.size());
  std::fill(lossmatrix.begin(), lossmatrix.end(), R_PosInf);
  int n = G0.rows();
  int p = G0.cols();
  Eigen::MatrixXd G(n, 2*p);

  double loss;
  double tau;
  double mu;

  if(dichotomous_i) {
    // no optimization needed
    loss = Gfitbeta_cpp(r, h, G0, x, param, var_epsilon, infeaturesfit, List::create(0,0), dichotomous_i, G);
    tau = 999.9;
    mu = 0;
  } else {
    for(int indexmu = 0; indexmu < mugridi.size(); indexmu++) {
      for(int indextau = 0; indextau < taugrid.size(); indextau++) {
        lossmatrix(indextau, indexmu) = Gfitbeta_cpp(r, h, G0, x, param, var_epsilon, infeaturesfit, List::create(mugridi[indexmu],log(taugrid[indextau])), dichotomous_i, G);
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

  result[0] = loss;
  result[1] = tau;
  result[2] = mu;

  return result;
}


// [[Rcpp::export]]
NumericMatrix loopfeatures_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0,
                               Eigen::MatrixXd x, NumericVector ifit, LogicalVector infeatures,
                               Eigen::MatrixXd mugrid, LogicalVector dichotomous, NumericVector taugrid,
                               List param, double var_epsilon) {
  int p = x.cols();
  // std::cout << "The value of x is: " << p << std::endl;
  NumericMatrix outputmatrix(p, 3);

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
    Eigen::VectorXd xi = x.col(i);
    LogicalVector infeaturesfit = updateinfeatures_cpp(infeatures, i);
    bool dichotomous_i = dichotomous[i];
    Eigen::VectorXd mugridi = mugrid.col(i);

    Eigen::VectorXd result = add_depth_cpp(xi, r, h, G0, infeaturesfit,
                                           dichotomous_i, mugridi, taugrid, param, var_epsilon);

    // assign the returned vector to a row in the output matrix
    outputmatrix(i, 0) = result[0];
    outputmatrix(i, 1) = result[1];
    outputmatrix(i, 2) = result[2];
  });

  return outputmatrix;

}















