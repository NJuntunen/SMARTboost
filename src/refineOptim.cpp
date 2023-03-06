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
#include <nlopt.hpp>
#include "sigmoidf_cpp.h"

using namespace Rcpp;

// Gfitbeta function
double Gfitbeta2_cpp(const Eigen::VectorXd r, const Eigen::VectorXd h, const Eigen::MatrixXd G0, const Eigen::VectorXd xi, SMARTParamStruct SMARTparams,
                    const double var_epsilon, const double muv, double tau,
                    const bool dichotomous_i) {

  const int n = G0.rows();
  const int p = G0.cols();
  Eigen::MatrixXd G(n, 2*p);

  Eigen::VectorXd gL(xi.size());

  double mu = muv;

  tau = std::max(tau, 0.2);


  gL = sigmoidf_cpp(xi, mu, tau, SMARTparams.sigmoid, dichotomous_i);

  G = updateG_allocated_cpp(G0, gL, G);

  FitBetaStruct result = fitbeta_cpp(r, G, var_epsilon, SMARTparams, mu, tau, dichotomous_i);

  double loss = result.loss;

  return loss;
}

double objective_function(const std::vector<double> &x, std::vector<double> &grad, void *data){

  RefineOptimData* obj_data = reinterpret_cast<RefineOptimData*>(data);

  Eigen::VectorXd &r = obj_data -> r;
  Eigen::VectorXd &h = obj_data -> h;
  Eigen::MatrixXd &G0 = obj_data -> G0;
  Eigen::VectorXd &xi = obj_data -> xi;
  SMARTParamStruct &SMARTparams = obj_data -> SMARTparams;
  double &var_epsilon = obj_data -> var_epsilon;
  double &tau = obj_data -> tau;
  bool &dichotomous_i = obj_data -> dichotomous_i;

  return Gfitbeta2_cpp(r, h, G0, xi, SMARTparams, var_epsilon, x[0],  tau, dichotomous_i);
}

std::vector<double> optimize_mutau_cpp (Eigen::VectorXd r,Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi,
                       SMARTParamStruct SMARTparams, double var_epsilon ,double tau, bool dichotomous_i,double mu0) {

  std::vector<double> res(2);
  const int n = G0.rows();
  const int p = G0.cols();
  Eigen::MatrixXd G(n, 2*p);
  double xtol_rel = SMARTparams.xtolOptim/(1+tau);

  RefineOptimData obj_data;
  obj_data.r = r;
  obj_data.h = h;
  obj_data.G0 = G0;
  obj_data.xi = xi;
  obj_data.SMARTparams = SMARTparams;
  obj_data.var_epsilon = var_epsilon;
  obj_data.dichotomous_i = dichotomous_i;
  obj_data.tau = tau;

  nlopt::opt opt(nlopt::LD_LBFGS, 1);
  opt.set_min_objective(objective_function, &obj_data);
  opt.set_maxeval(100);
  opt.set_xtol_rel(xtol_rel);
  std::vector<double> x(1);
  x[0] = mu0;
  double minf;

  try{
    nlopt::result result = opt.optimize(x, minf);
  }catch(std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }

  res[0] = minf;
  res[1] = x[0];

  return res;
}

std::vector<double>  refineOptim_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi,
                                Eigen::VectorXd dichotomous, double mu0, bool dichotomous_i,double tau0, SMARTParamStruct SMARTparams, double var_epsilon) {

  std::vector<double> result(3);
  std::vector<double> taugrid;
  Eigen::MatrixXd output(7,2);
  double loss;
  double tau;
  double mu;

  if (dichotomous_i) {
    Eigen::VectorXd gL = sigmoidf_cpp(xi,mu0,tau0,SMARTparams.sigmoid, dichotomous_i);
    const int n = G0.rows();
    const int p = G0.cols();
    Eigen::MatrixXd G(n, 2*p);
    G = updateG_allocated_cpp(G0,gL,G);
    double loss = std::numeric_limits<double>::infinity();
    double tau = tau0;
    double mu = mu0;
  } else {
    if (SMARTparams.sharptree) {

      taugrid.push_back(tau0);

    } else if (SMARTparams.taugridpoints == 1) {

      double start = -2.7;
      double end = 2.7;
      double step = 0.3;

      for (double i = start; i <= end; i += step) {
        double value = std::exp(std::log(tau0) + i);
        taugrid.push_back(value);
      }
    } else if (SMARTparams.taugridpoints == 2) {

      double start = -1.8;
      double end = 1.8;
      double step = 0.3;

      for (double i = start; i <= end; i += step) {
        double value = std::exp(std::log(tau0) + i);
        taugrid.push_back(value);
      }

    } else {
      if (tau0 < 8.0) {

        double start = -0.9;
        double end = 0.9;
        double step = 0.3;

        for (double i = start; i <= end; i += step) {
          double value = std::exp(std::log(tau0) + i);
          taugrid.push_back(value);
        }
      } else {

        double start = -0.9;
        double end = 1.8;
        double step = 0.3;

        for (double i = start; i <= end; i += step) {
          double value = std::exp(std::log(tau0) + i);
          taugrid.push_back(value);
        }
      }
    }

    int loops = taugrid.size();
    std::vector<double> res(2);
    std::mutex mtx;
    RcppThread::parallelFor(0, taugrid.size(), [&](int i) {
    // for (int j = 0; j < loops; j++) {

      double tau = taugrid[i];

      res = optimize_mutau_cpp(r,h,G0,xi,SMARTparams,var_epsilon,tau,dichotomous_i,mu0);

      mtx.lock();
      output(i,0) = res[0];
      output(i,1) = res[1];
      mtx.unlock();
    });

  }

  int minindex;
  double minValue = output.col(0).minCoeff(&minindex);

  loss = output(minindex,0);
  tau = taugrid[minindex];
  mu = output(minindex,1);

  result[0] = loss;
  result[1] = tau;
  result[2] = mu;

  return result;
}



