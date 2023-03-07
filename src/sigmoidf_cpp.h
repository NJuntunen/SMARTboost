#include <Rcpp.h>
#include "structures.h"

Eigen::VectorXd sigmoidf_cpp(Eigen::VectorXd x, double mu, double tau, std::string sigmoid, bool dichotomous);
Eigen::MatrixXd updateG_allocated_cpp(Eigen::MatrixXd G0, Eigen::VectorXd g, Eigen::MatrixXd G);
FitBetaStruct fitbeta_cpp(Eigen::VectorXd r, Eigen::MatrixXd G, double var_epsilon, SMARTParamStruct SMARTparams, double mu, double tau, bool dichotomous_i);
Eigen::MatrixXd loopfeatures_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0,
                                 Eigen::MatrixXd x, Eigen::MatrixXd mugrid, Eigen::VectorXd dichotomous, Eigen::VectorXd taugrid,
                                 SMARTParamStruct SMARTparams, const double var_epsilon);
std::vector<double>  refineOptim_cpp(Eigen::VectorXd r, Eigen::VectorXd h, Eigen::MatrixXd G0, Eigen::VectorXd xi,
                                     Eigen::VectorXd dichotomous, double mu0, bool dichotomous_i,double tau0,SMARTParamStruct SMARTparams, double var_epsilon);
double objective_function(const std::vector<double> &x, std::vector<double> &grad, void *data);
double Gfitbeta2_cpp(const Eigen::VectorXd r, const Eigen::VectorXd h, const Eigen::MatrixXd G0, const Eigen::VectorXd xi, SMARTParamStruct SMARTparams,
                     const double var_epsilon, const double muv, double tau,
                     const bool dichotomous_i);
