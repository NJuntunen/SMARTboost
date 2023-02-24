#include <Rcpp.h>
#include <nlopt.hpp>
#include <cmath>
#include <math.h>

// Define the objective function
double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  double a = x[0], b = x[1];
  if (!grad.empty()) {
    grad[0] = 2 * a;
    grad[1] = 2 * b;
  }
  return a*a + b*b;
}

// [[Rcpp::export]]
Rcpp::List nlopt_example() {
  nlopt::opt opt(nlopt::LD_MMA, 2);
  std::vector<double> lb(2);
  lb[0] = -std::numeric_limits<double>::infinity(); lb[1] = -std::numeric_limits<double>::infinity();
  opt.set_lower_bounds(lb);
  opt.set_min_objective(myfunc, NULL);
  opt.set_xtol_rel(1e-4);
  std::vector<double> x(2);
  x[0] = 1.234; x[1] = 5.678;
  double minf;
  nlopt::result result = opt.optimize(x, minf);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["optimal_value"] = minf,
    Rcpp::_["optimal_solution"] = Rcpp::NumericVector(x.begin(), x.end()),
    Rcpp::_["result_code"] = x
  );
  return out;
}
