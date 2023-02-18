// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <RcppThread.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// Define the struct for the parameters
struct ParamStruct {
  String loss;
  double coeff;
  String verbose;
  bool randomizecv;
  int ncores;
  double sharevalidation;
  double stderulestop;
  bool stopwhenlossup;
  double lambda;
  int depth;
  String sigmoid;
  double meanlntau;
  double varlntau;
  double doflntau;
  double varmu;
  double dofmu;
  double subsamplesharevs;
  bool subsamplefinalbeta;
  double subsampleshare_columns;
  double mugridpoints;
  double taugridpoints;
  bool refineOptimGrid;
  double xtolOptim;
  bool optimizevs;
  bool sharptree;
  int ntrees;
  double R2p;
  double p0;
  double loglikdivide;
  double overlap;
};

// // sigmoidf function
// // [[Rcpp::export]]
// Eigen::VectorXd sigmoidf_cpp(Eigen::VectorXd x, double mu, double tau, std::string sigmoid, bool dichotomous) {
//
//   Eigen::VectorXd g(x.size());
//   if(dichotomous){
//     for (int i = 0; i < x.size(); i++) {
//       g(i) = (x(i) > 0) ? 1 : 0;
//     }
//     return g;
//   } else {
//     if(sigmoid == "sigmoidsqrt"){
//       g = 0.5 + 0.5 * ( 0.5 * tau * (x.array() - mu) / sqrt((1 + pow(0.5 * tau * (x.array() - mu),2))));
//     } else if(sigmoid == "sigmoidlogistic"){
//       g = 1 - 1/(1 + exp(tau * (x.array() - mu)));
//     }
//   }
//
//   return g;
// }
//
// // updateG_allocated function
// // [[Rcpp::export]]
// Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, Eigen::VectorXd g, Eigen::MatrixXd G) {
//   int n = G0.rows(), p = G0.cols();
//   for (int i = 0; i < p; i++) {
//     for (int j = 0; j < n; j++) {
//       G(j, i) = G0(j, i) * g[j];
//       G(j, i + p) = G0(j, i) * (1 - g[j]);
//     }
//   }
//   return G;
// }
//
// // define lnpmu C++ function
// double lnpmu(double mu, double varmu, double dofmu) {
//   double s = sqrt(varmu);
//   double lnp = R::dnorm(mu, 0.0, s, true);
//   return lnp;
// }
//
// // logpdft function
// // x: input values
// // m: mean
// // s: standard deviation
// // v: degrees of freedom
// NumericVector logpdft(NumericVector x, double m, double s, double v) {
//   double constant = -0.5723649429247001 + lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v);
//   NumericVector z = (x - m) / s;
//   NumericVector logpdfz = constant - 0.5 * (1 + v) * log(1 + pow(z, 2) / v);
//   return logpdfz - log(s);
// }
//
// // lnptau function
// // tau: input values
// // meanlntau: mean of log(tau)
// // varlntau: variance of log(tau)
// // doflntau: degrees of freedom
// // depth: tree depth
// double lnptau(double tau, double meanlntau, double varlntau, double doflntau, int depth) {
//   double s = sqrt(varlntau / depth);
//   NumericVector logpdf_tau = logpdft(log(tau), meanlntau, s, doflntau);
//   double lnp_tau = sum(logpdf_tau);
//   return lnp_tau;
// }
//
//
// // define the main function
// // [[Rcpp::export]]
// List fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, List param, LogicalVector infeaturesfit,
//                  double mu, double tau, bool dichotomous_i) {
//
//   Eigen::VectorXd diagGGh = G.colwise().squaredNorm();
//   Eigen::MatrixXd GGh = G.transpose() * G;
//   int n = G.rows();
//   int p = G.cols();
//   double var_r = var_epsilon/(1-as<double>(param["R2p"]));
//   double Pb = diagGGh.sum()/(n*var_r*as<double>(param["R2p"]));
//   Eigen::MatrixXd I_p = Eigen::MatrixXd::Identity(p,p);
//   Eigen::MatrixXd GGh_var_r_Pb_I_p = GGh + var_r*as<double>(param["loglikdivide"])*Pb*I_p;
//   Eigen::VectorXd beta(p);
//   List ret;
//
//   // // Use BiCGSTAB to solve the system
//   // Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
//   // solver.compute(*GGh_var_r_Pb_I_p);
//   // *beta = solver.solve(G.transpose() * r);
//   //
//   // // Check if the solution is valid, and if not, use a larger Pb
//   // int max_iter = 5; // set a maximum number of iterations
//   // int iter = 0;
//   // while (solver.info() != Eigen::Success && iter < max_iter) {
//   //   iter++;
//   //   std::unique_ptr<double> new_Pb(new double((*Pb) * 2.01)); // create a new unique_ptr to avoid memory leaks
//   //   *GGh_var_r_Pb_I_p = *GGh + var_r * as<double>(param["loglikdivide"]) * (*new_Pb) * (*I_p);
//   //   solver.compute(*GGh_var_r_Pb_I_p);
//   //   *beta = solver.solve(G.transpose() * r);
//   //   std::swap(Pb, new_Pb); // release old memory and assign new value to Pb
//   // }
//   // if (solver.info() != Eigen::Success) {
//   //   Rcout << "Warning: solver did not converge\n";
//   // }
//   //
//   // std::unique_ptr<Eigen::MatrixXd> Gbeta(new Eigen::MatrixXd(G * (*beta)));
//   // std::unique_ptr<double> loglik(new double(-0.5*(r-(*Gbeta)).squaredNorm()/var_r/as<double>(param["loglikdivide"])));
//   // std::unique_ptr<double> logpdfbeta(new double(-0.5*(p*log(2*M_PI) - p*log((*Pb)) + *Pb*(*beta).squaredNorm())));
//   //
//   // double logpdfmu = 0;
//   // double logpdftau = 0;
//   //
//   // if (!dichotomous_i) {
//   //   if (as<double>(param["sharptree"])) {
//   //     logpdfmu = lnpmu(mu, as<double>(param["varmu"]), as<double>(param["dofmu"]));
//   //   } else {
//   //     logpdfmu = lnpmu(mu, as<double>(param["varmu"]), as<double>(param["dofmu"]));
//   //     logpdftau = lnptau(tau, as<double>(param["meanlntau"]), as<double>(param["varlntau"]), as<double>(param["doflntau"]), as<double>(param["depth"]));
//   //   }
//   // }
//   //
//   // double loss = -(*loglik + (*logpdfbeta) + logpdftau + logpdfmu);
//
//   ret["loss"] = 10;
//   ret["Gbeta"] = 10;
//   ret["beta"] = 10;
//
//   return ret;
// }
//
//
//
// // Gfitbeta function
// double Gfitbeta_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi, List param,
//                     double var_epsilon, LogicalVector infeaturesfit, List muLogTau,
//                     bool dichotomous_i, Eigen::MatrixXd G) {
//   Eigen::VectorXd gL(xi.size());
//   double tau = exp(as<double>(muLogTau[1]));
//   tau = std::max(tau, 0.2);
//
//
//   // gL = sigmoidf_cpp(xi, as<double>(muLogTau[0]), tau, param["sigmoid"], dichotomous_i);
//
//   // G = updateG_allocated_cpp(G0, gL, G);
//
//   // List result = fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit, as<double>(muLogTau[0]), tau, dichotomous_i);
//
//
//   return 10;
// }

RcppParallel::RVector<int> updateinfeatures_cpp(RcppParallel::RVector<int> infeatures, int ifit) {

  infeatures[ifit] = 1;

  return infeatures;
}

// // Function that calculates the loss, tau, and mu for one column of x
// std::vector<double> add_depth_cpp(RcppParallel::RMatrix<double>::Column x, RcppParallel::RVector<double> r, RcppParallel::RVector<double> h,
//                                             RcppParallel::RMatrix<double> G0, RcppParallel::RVector<int> infeaturesfit,
//                                             bool dichotomous_i, RcppParallel::RMatrix<double>::Column mugridi,
//                               RcppParallel::RVector<double> taugrid, ParamStruct SMARTparams, double var_epsilon) {
//
//   std::vector<double> result(3, 0);
//
//
//   const int nrow = taugrid.length();
//   const int ncol = mugridi.size();
//
//   std::vector<std::vector<double>> lossmatrix(nrow, std::vector<double>(ncol, R_PosInf));  // create a 2D vector of size nrow x ncol with default value R_PosInf
//
//   int n = G0.nrow();
//   int p = G0.ncol();
//
//   double loss;
//   double tau;
//   double mu;
//
//   // if(dichotomous_i) {
//   //   // no optimization needed
//   //   loss = Gfitbeta_cpp(r, h, G0, x, param, var_epsilon, infeaturesfit, List::create(0,0), dichotomous_i, G);
//   //   tau = 999.9;
//   //   mu = 0;
//   // } else {
//   //   for(int indexmu = 0; indexmu < mugridi.size(); indexmu++) {
//   //     for(int indextau = 0; indextau < taugrid.size(); indextau++) {
//   //       lossmatrix(indextau, indexmu) = Gfitbeta_cpp(r, h, G0, x, param, var_epsilon, infeaturesfit, List::create(mugridi[indexmu],log(taugrid[indextau])), dichotomous_i, G);
//   //       if (indextau > 0 && lossmatrix(indextau, indexmu) > lossmatrix(indextau-1, indexmu)) {
//   //         break;
//   //       }
//   //     }
//   //   }
//   //
//   //   int min_row = 0;
//   //   int min_col = 0;
//   //   double min_val = lossmatrix(0, 0);
//   //
//   //   for (int i = 0; i < lossmatrix.nrow(); i++) {
//   //     for (int j = 0; j < lossmatrix.ncol(); j++) {
//   //       double val = lossmatrix(i, j);
//   //       if (val < min_val) {
//   //         min_val = val;
//   //         min_row = i;
//   //         min_col = j;
//   //       }
//   //     }
//   //   }
//   //
//   //   loss = lossmatrix(min_row, min_col);
//   //   tau = taugrid(min_row);
//   //   mu = mugridi(min_col);
//   //
//   //   // Optionally, further optimize over mu
//   //   if(param["optimizevs"]) {
//   //
//   //     //define nlopt optimization here
//   //
//   //   }
//   // }
//   loss = 10;
//   tau = 10;
//   mu = 10;
//   result[0] = loss;
//   result[1] = tau;
//   result[2] = mu;
//
//   return result;
// }
//

// [[Rcpp::export]]
NumericMatrix loopfeatures_cpp(NumericVector r, NumericVector h, NumericMatrix G0,
                               NumericMatrix x, IntegerVector ifit, IntegerVector infeatures,
                               NumericMatrix mugrid, LogicalVector dichotomous, NumericVector taugrid,
                               List param, double var_epsilon) {

  //Defining RcppParallel objects for parallel loop
  RcppParallel::RVector<double> rp(r);
  RcppParallel::RVector<double> hp(h);
  RcppParallel::RMatrix<double> G0p(G0);
  RcppParallel::RMatrix<double> xp(x);
  RcppParallel::RVector<int> ifitp(ifit);
  RcppParallel::RVector<int> infeaturesp(infeatures);
  RcppParallel::RMatrix<double> mugridp(mugrid);
  RcppParallel::RVector<int> dichotomousp(dichotomous);
  RcppParallel::RVector<double> taugridp(taugrid);

//Defining parameters
ParamStruct SMARTparams;
SMARTparams.loss = as<String>(param["loss"]);
SMARTparams.coeff = as<double>(param["coeff"]);
SMARTparams.verbose = as<String>(param["verbose"]);
SMARTparams.randomizecv = as<bool>(param["randomizecv"]);
SMARTparams.ncores = as<int>(param["ncores"]);
SMARTparams.stderulestop = as<double>(param["stderulestop"]);
SMARTparams.lambda = as<double>(param["lambda"]);
SMARTparams.sigmoid = as<String>(param["sigmoid"]);
SMARTparams.meanlntau = as<double>(param["meanlntau"]);
SMARTparams.varlntau = as<double>(param["varlntau"]);
SMARTparams.doflntau = as<double>(param["doflntau"]);
SMARTparams.varmu = as<double>(param["varmu"]);
SMARTparams.dofmu = as<double>(param["dofmu"]);
SMARTparams.subsamplesharevs = as<double>(param["subsamplesharevs"]);
SMARTparams.subsamplefinalbeta = as<bool>(param["subsamplefinalbeta"]);
SMARTparams.subsampleshare_columns = as<double>(param["subsampleshare_columns"]);
SMARTparams.refineOptimGrid = as<bool>(param["refineOptimGrid"]);
SMARTparams.xtolOptim = as<double>(param["xtolOptim"]);
SMARTparams.optimizevs = as<bool>(param["optimizevs"]);
SMARTparams.sharptree = as<bool>(param["sharptree"]);
SMARTparams.R2p = as<double>(param["R2p"]);
SMARTparams.p0 = as<double>(param["p0"]);
SMARTparams.loglikdivide = as<double>(param["loglikdivide"]);
SMARTparams.overlap = as<double>(param["overlap"]);

  int p = x.cols();
//   // std::cout << "The value of x is: " << p << std::endl;
  NumericMatrix outputmatrix(p, 3);
  RcppParallel::RMatrix<double> output(outputmatrix);
  IntegerVector ps(p);
  RcppParallel::RVector<int> psp(ps);
  for (int i = 0; i < p; i++) {
    psp[i] = i+1;
  }

  if (SMARTparams.subsampleshare_columns < 1) {
    int psmall = round(p * SMARTparams.subsampleshare_columns);
    IntegerVector p_sampled = Rcpp::sample(ps, psmall);
    ps = p_sampled;
  }

  // loop through the columns in parallel
  std::vector<double> xi(x.nrow());
  RcppThread::parallelFor(0, p, [&](int i) {

    int index = ps[i]-1;
    RcppParallel::RMatrix<double>::Column xip = xp.column(index);
    RcppParallel::RVector<int> infeaturesfit = updateinfeatures_cpp(infeaturesp, index);
    bool dichotomous_i = dichotomous[index];
    RcppParallel::RMatrix<double>::Column mugridip = mugridp.column(index);
//
//     std::vector<double> result = add_depth_cpp(xip, rp, hp, G0p, infeaturesfit,
//                                                dichotomous_i, mugridip, taugridp, SMARTparams, var_epsilon);

    // assign the returned vector to a row in the output matrix
    output(i, 0) = 10;
    output(i, 1) = 10;
    output(i, 2) = 10;
  });

  // Create a new NumericMatrix to hold the data
  NumericMatrix Routput(output.nrow(), output.ncol());

  // Copy the data from the RMatrix to the NumericMatrix
  for(int i = 0; i < output.nrow(); i++) {
    for(int j = 0; j < output.ncol(); j++) {
      Routput(i, j) = Routput(i, j);
    }
  }

  return Routput;

}











