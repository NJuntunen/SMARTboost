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

        Eigen::VectorXi ssi_indices = Eigen::Map<Eigen::VectorXi>(ssi.begin(), ssi.size());

        Eigen::VectorXd r_subsampled = r.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::VectorXd h_subsampled = h.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::MatrixXd G0_subsampled = G0.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, G0.cols());
        Eigen::MatrixXd x_subsampled = x.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, x.cols());

        outputarray = loopfeatures_cpp(r_subsampled, h_subsampled, G0_subsampled, x_subsampled, mugrid, dichotomous, taugrid, SMARTparams, var_epsilon); // loops over all variables
      }else {

        Eigen::VectorXi ssi_indices = Eigen::Map<Eigen::VectorXi>(ssi.begin(), ssi.size());
        Eigen::VectorXd r_subsampled = r.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::MatrixXd G0_subsampled = G0.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, G0.cols());
        Eigen::MatrixXd x_subsampled = x.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, x.cols());

        outputarray = loopfeatures_cpp(r_subsampled, h, G0_subsampled, x_subsampled, mugrid, dichotomous, taugrid, SMARTparams, var_epsilon);
      }
    }

    int i = (outputarray.col(0).array() == outputarray.col(0).minCoeff()).cast<int>().maxCoeff();
    double tau0 = outputarray(i, 1);
    double mu0 = outputarray(i, 2);
    std::vector<double> opt_result(3);

    std::cout << "start optim " << std::endl;
    if (SMARTparams.subsamplesharevs < 1 && SMARTparams.subsamplefinalbeta == true) {
      if (h.size() == 1) {
        std::cout << "opt 1" << std::endl;
        Eigen::VectorXi ssi_indices = Eigen::Map<Eigen::VectorXi>(ssi.begin(), ssi.size());
        Eigen::VectorXd r_subsampled = r.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::VectorXd h_subsampled = h.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::MatrixXd G0_subsampled = G0.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, G0.cols());
        Eigen::MatrixXd x_subsampled = x.block(ssi_indices.minCoeff() - 1, i, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, 1);

        opt_result = refineOptim_cpp(r_subsampled, h_subsampled, G0_subsampled, x_subsampled, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);
      } else {
        std::cout << "opt 2" << std::endl;
        Eigen::VectorXi ssi_indices = Eigen::Map<Eigen::VectorXi>(ssi.begin(), ssi.size());
        Eigen::VectorXd r_subsampled = r.segment(ssi_indices.minCoeff() - 1, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1);
        Eigen::MatrixXd G0_subsampled = G0.block(ssi_indices.minCoeff() - 1, 0, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, G0.cols());
        Eigen::MatrixXd x_subsampled = x.block(ssi_indices.minCoeff() - 1, i, ssi_indices.maxCoeff() - ssi_indices.minCoeff() + 1, 1);

        opt_result = refineOptim_cpp(r_subsampled, h, G0_subsampled, x_subsampled, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);
      }
    } else {

      Eigen::VectorXd xi = x.col(i);

      opt_result = refineOptim_cpp(r, h, G0, xi, dichotomous, mu0, dichotomous[i], tau0, SMARTparams, var_epsilon);

    }
    std::cout << "optim done " << std::endl;
    double loss = opt_result[0];
    double tau = opt_result[1];
    double mu = opt_result[2];

    Eigen::VectorXd gL;
    Eigen::MatrixXd G(n, 2^(depth+1));
    gL = sigmoidf_cpp(x.col(i), mu, tau, SMARTparams.sigmoid, dichotomous[i]);
    std::cout << "sigmoid done " << std::endl;

    try {
      // code that may potentially cause an error
      G = updateG_allocated_cpp(G0, gL, G);
    } catch (std::exception& e) {
      // handle the error and return an error message
      Rcout << "Error: " << e.what() << std::endl;
      break;
    }


    std::cout << "G done " << std::endl;

    FitBetaStruct FitBeta;
    Eigen::VectorXd yfit;
    Eigen::VectorXd beta;

    FitBeta = fitbeta_cpp(r, G, var_epsilon, SMARTparams, mu, tau, dichotomous[i]);
    yfit = FitBeta.Gbeta;
    beta = FitBeta.beta;
    loss = FitBeta.loss;
    std::cout << "fitbeta done " << std::endl;

    fi2[depth] = (pow(yfit.array(),2).sum() - pow(yfit0.array(),2).sum()) / n;

    G0 = G;
    loss0 = loss;
    yfit0 = yfit;
    ifit.conservativeResize(1);
    std::cout << "ifit size: " << ifit.size() << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(5));
    ifit(ifit.size()-1) = i;
    mufit.conservativeResize(1);
    mufit(mufit.size()-1) = mu;
    taufit.conservativeResize(1);
    taufit(taufit.size()-1) = tau;
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
