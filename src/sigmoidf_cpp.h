#include <Rcpp.h>
#include "structures.h"
Eigen::VectorXd sigmoidf_cpp(Eigen::VectorXd x, double mu, double tau, std::string sigmoid, bool dichotomous);
Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, Eigen::VectorXd g, Eigen::MatrixXd G);
FitBetaStruct fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, SMARTParamStruct SMARTparams, double mu, double tau, bool dichotomous_i);
