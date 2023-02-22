#include <RcppEigen.h>
#include <Rcpp.h>
// Define the struct for the parameters
struct SMARTParamStruct {
  std::string loss;
  double coeff;
  std::string verbose;
  bool randomizecv;
  int ncores;
  double sharevalidation;
  double stderulestop;
  bool stopwhenlossup;
  double lambda;
  int depth;
  std::string sigmoid;
  double meanlntau;
  double varlntau;
  double doflntau;
  double varmu;
  double dofmu;
  double subsamplesharevs;
  bool subsamplefinalbeta;
  double subsampleshare_columns;
  double mugridpoints;
  int taugridpoints;
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

struct FitBetaStruct {
  double loss;
  Eigen::VectorXd Gbeta;
  Eigen::VectorXd beta;

};

struct RefineOptimData {
  Eigen::VectorXd r;
  Eigen::VectorXd h;
  Eigen::MatrixXd G0;
  Eigen::VectorXd xi;
  SMARTParamStruct SMARTparams;
  double var_epsilon;
  double tau;
  bool dichotomous_i;
  Eigen::MatrixXd G;

};



