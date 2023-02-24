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
                    const double var_epsilon, const std::vector<double> muv, double tau,
                    const bool dichotomous_i, Eigen::MatrixXd G) {

  Eigen::VectorXd gL(xi.size());

  double mu = muv[1];

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
  Eigen::MatrixXd &G = obj_data -> G;

  double result = Gfitbeta2_cpp(r, h, G0, xi, SMARTparams, var_epsilon, x,  tau, dichotomous_i, G);

  return result;
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

  nlopt::opt opt(nlopt::LD_MMA,1);
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
  std::cout << "opt done" << std::endl;

  res[0] = minf;
  res[1] = x[0];
  // double result = Gfitbeta2_cpp(r, h, G0, xi, SMARTparams, var_epsilon, x,  tau, dichotomous_i, G);
  // res[0] = result;
  // res[1] = 8;


  return res;
}


// [[Rcpp::export]]
Eigen::MatrixXd refineOptim_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi,
                                Eigen::VectorXd dichotomous, double mu0, bool dichotomous_i,double tau0,List param, double var_epsilon) {

  SMARTParamStruct SMARTparams;
  SMARTparams.randomizecv = as<bool>(param["randomizecv"]);
  SMARTparams.ncores = as<int>(param["ncores"]);
  SMARTparams.lambda = as<double>(param["lambda"]);
  SMARTparams.depth = as<int>(param["depth"]);
  SMARTparams.sigmoid = as<String>(param["sigmoid"]);
  SMARTparams.meanlntau = as<double>(param["meanlntau"]);
  SMARTparams.varlntau = as<double>(param["varlntau"]);
  SMARTparams.doflntau = as<double>(param["doflntau"]);
  SMARTparams.varmu = as<double>(param["varmu"]);
  SMARTparams.dofmu = as<double>(param["dofmu"]);
  SMARTparams.subsamplesharevs = as<double>(param["subsamplesharevs"]);
  SMARTparams.subsampleshare_columns = as<double>(param["subsampleshare_columns"]);
  SMARTparams.xtolOptim = as<double>(param["xtolOptim"]);
  SMARTparams.optimizevs = as<bool>(param["optimizevs"]);
  SMARTparams.sharptree = as<bool>(param["sharptree"]);
  SMARTparams.R2p = as<double>(param["R2p"]);
  SMARTparams.p0 = as<double>(param["p0"]);
  SMARTparams.loglikdivide = as<double>(param["loglikdivide"]);
  SMARTparams.overlap = as<double>(param["overlap"]);
  SMARTparams.taugridpoints = as<int>(param["taugridpoints"]);

  std::vector<double> result(3);
  std::vector<double> taugrid;
  Eigen::MatrixXd output(7,2);


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
    std::cout << taugrid.size() << std::endl;
    int loops = taugrid.size();
    std::vector<double> res(2);
    //RcppThread::parallelFor(0, taugrid.size(), [&](int i) {
    for (int j = 0; j < loops; j++) {

      double tau = taugrid[j];

      res = optimize_mutau_cpp(r,h,G0,xi,SMARTparams,var_epsilon,tau,dichotomous_i,mu0);

      output(j,0) = res[0];
      output(j,1) = res[1];

      // output(j,0) = 10;
      // output(j,1) = 10;
      std::cout << j << std::endl;

    }


  }
  //
  //         minindex <- which.min(lossmatrix[,1])
  //         loss <- lossmatrix[minindex,1]
  //       tau <- taugrid[minindex]
  //       mu <- lossmatrix[minindex,2]

  return output;
}



