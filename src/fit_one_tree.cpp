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

// [[Rcpp::export]]
List fit_one_tree_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd x, Eigen::MatrixXd mugrid,
                      Eigen::VectorXd dichotomous, Eigen::VectorXd taugrid, List param) {


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
  SMARTparams.subsamplefinalbeta = as<bool>(param["subsamplefinalbeta"]);
  SMARTparams.xtolOptim = as<double>(param["xtolOptim"]);
  SMARTparams.optimizevs = as<bool>(param["optimizevs"]);
  SMARTparams.sharptree = as<bool>(param["sharptree"]);
  SMARTparams.R2p = as<double>(param["R2p"]);
  SMARTparams.p0 = as<double>(param["p0"]);
  SMARTparams.loglikdivide = as<double>(param["loglikdivide"]);
  SMARTparams.overlap = as<double>(param["overlap"]);
  SMARTparams.taugridpoints = as<int>(param["taugridpoints"]);

  double mean = r.mean();
  Eigen::VectorXd deviations = r.array() - mean;
  double squaredDeviations = deviations.array().square().sum();
  double var_wr = squaredDeviations / r.size();

  double var_epsilon = var_wr * (1- SMARTparams.R2p);

  int n = x.rows();
  int p = x.cols();

  Eigen::MatrixXd G0(n, 1);
  G0.setOnes();
  double loss0 = std::numeric_limits<double>::infinity();
  Eigen::VectorXd ifit;
  Eigen::VectorXd mufit;
  Eigen::VectorXd taufit;
  Eigen::VectorXd betafit;
  Eigen::VectorXd yfit0(n);
  yfit0.setZero();
  Eigen::VectorXd fi2(SMARTparams.depth);
  fi2.setZero();

  double subsamplesize = round(n * SMARTparams.subsamplesharevs);

  IntegerVector ssi(n);

  if (SMARTparams.subsamplesharevs == 1) {

    for (int i = 0; i < n; i++) {
      ssi[i] = i;
    }
  }else {
    IntegerVector n_sampled = Rcpp::sample(n, subsamplesize);
    ssi = n_sampled;
  }

  for(int depth = 0; depth < SMARTparams.depth; depth++){

    int cols = p;
    if (SMARTparams.subsampleshare_columns < 1) {
      int psmall = round(p * SMARTparams.subsampleshare_columns);
      cols = psmall;
    }

    Eigen::MatrixXd outputarray(cols, 3);

    if(SMARTparams.subsamplesharevs == 1){

      outputarray = loopfeatures_cpp(r, h, G0, x, mugrid, dichotomous, taugrid, SMARTparams, var_epsilon);

    }else { //Variable selection using a random sub-set of the sample. All the sample is then used in refinement.
      if (h.size() == 1){

        Eigen::VectorXd r_subsampled(ssi.size());
        Eigen::VectorXd h_subsampled(ssi.size());
        Eigen::MatrixXd G0_subsampled(ssi.size(), G0.cols());
        Eigen::MatrixXd x_subsampled(ssi.size(), x.cols());
        int count = 0;
        for(int i : ssi){

          r_subsampled(count) = r(i);
          h_subsampled(count) = h(i);
          G0_subsampled.row(count) = G0.row(count);
          x_subsampled.row(count) = x.row(count);
          count = count + 1;
        }

        outputarray = loopfeatures_cpp(r_subsampled, h_subsampled, G0_subsampled, x_subsampled, mugrid, dichotomous, taugrid, SMARTparams, var_epsilon); // loops over all variables
      }else {

        Eigen::VectorXd r_subsampled(ssi.size());
        Eigen::MatrixXd G0_subsampled(ssi.size(), G0.cols());
        Eigen::MatrixXd x_subsampled(ssi.size(), x.cols());
        int count = 0;
        for(int i : ssi){

          r_subsampled(count) = r(i);
          G0_subsampled.row(count) = G0.row(count);
          x_subsampled.row(count) = x.row(count);
          count = count + 1;
        }

        outputarray = loopfeatures_cpp(r_subsampled, h, G0_subsampled, x_subsampled, mugrid, dichotomous, taugrid, SMARTparams, var_epsilon);
      }
    }

    int i = (outputarray.col(0).array() == outputarray.col(0).minCoeff()).cast<int>().maxCoeff();
    double tau0 = outputarray(i, 1);
    double mu0 = outputarray(i, 2);
    std::vector<double> opt_result(3);

    if (SMARTparams.subsamplesharevs < 1 && SMARTparams.subsamplefinalbeta == true) {
      if (h.size() == 1) {

        Eigen::VectorXd r_subsampled(ssi.size());
        Eigen::VectorXd h_subsampled(ssi.size());
        Eigen::MatrixXd G0_subsampled(ssi.size(), G0.cols());
        Eigen::VectorXd x_subsampled(ssi.size());
        int count = 0;
        for(int i : ssi){

          r_subsampled(count) = r(i);
          h_subsampled(count) = h(i);
          G0_subsampled.row(count) = G0.row(count);
          x_subsampled(count) = x(count,i);
          count = count + 1;
        }

        opt_result = refineOptim_cpp(r_subsampled, h_subsampled, G0_subsampled, x_subsampled, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);

      } else {

        Eigen::VectorXd r_subsampled(ssi.size());
        Eigen::MatrixXd G0_subsampled(ssi.size(), G0.cols());
        Eigen::VectorXd x_subsampled(ssi.size());
        int count = 0;
        for(int i : ssi){

          r_subsampled(count) = r(i);
          G0_subsampled.row(count) = G0.row(count);
          x_subsampled(count) = x(count,i);
          count = count + 1;
        }

        opt_result = refineOptim_cpp(r_subsampled, h, G0_subsampled, x_subsampled, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);
      }
    } else {

      Eigen::VectorXd xi = x.col(i);

      opt_result = refineOptim_cpp(r, h, G0, xi, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);

    }

    double loss = opt_result[0];
    double tau = opt_result[1];
    double mu = opt_result[2];

    int result = static_cast<int>(std::pow(static_cast<double>(2), static_cast<double>(depth + 1)));

    Eigen::VectorXd gL(x.rows());
    Eigen::MatrixXd G(n, result);

    gL = sigmoidf_cpp(x.col(i), mu, tau, SMARTparams.sigmoid, dichotomous[i]);
    G = updateG_allocated_cpp(G0, gL, G);

    FitBetaStruct FitBeta;
    Eigen::VectorXd yfit;
    Eigen::VectorXd beta;

    FitBeta = fitbeta_cpp(r, G, var_epsilon, SMARTparams, mu, tau, dichotomous[i]);
    yfit = FitBeta.Gbeta;
    beta = FitBeta.beta;
    loss = FitBeta.loss;

    fi2[depth] = (pow(yfit.array(),2).sum() - pow(yfit0.array(),2).sum()) / n;

    G0.resize(n, G.cols());
    G0 = G;
    loss0 = loss;
    yfit0 = yfit;
    ifit.conservativeResize(ifit.size() + 1);
    ifit(ifit.size() - 1) = i;
    std::cout << ifit << std::endl;
    mufit.conservativeResize(mufit.size() + 1);
    mufit(mufit.size() - 1) = mu;
    taufit.conservativeResize(taufit.size() + 1);
    taufit(taufit.size() - 1) = tau;
    betafit = beta;

  }

  Rcpp::List retlist;

  retlist["yfit0"] = yfit0;
  retlist["ifit"] = ifit;
  retlist["mufit"] = mufit;
  retlist["taufit"] = taufit;
  retlist["betafit"] = betafit;
  retlist["fi2"] = fi2;


  return retlist;
}
