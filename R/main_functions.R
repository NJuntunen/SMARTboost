#' @export
SMARTboost_fit <- function(x, ...) {
  UseMethod("SMARTboost_fit")
}

#' @export
#' @rdname SMARTboost_fit
SMARTboost_fit.recipe <- function(x, data, SMARTboost_model = NULL, param = SMARTparam(), ...) {

  processed <- hardhat::mold(x, data)

  default_param <- SMARTparam()
  new_param <- do.call(SMARTparam, list(...))
  new_param <- new_param[
    mapply(
      function(x, y) ifelse(is.null(x), !is.null(y), x != y),
      default_param,
      new_param)
  ]
  param <- utils::modifyList(param, as.list(new_param))

  SMARTboost_bridge(processed, param = param, SMARTboost_model)
}




SMARTboost_bridge <- function(processed, param = SMARTparam(), SMARTboost_model) {

  predictors <- processed$predictors
  outcomes <- processed$outcomes

  # initialize SMARTtrees
  data <- preparedataSMART(predictors, param)
  meanx <- data$meanx
  stdx <- data$stdx
  data <- data$data_standardized

  grids <- preparegridsSMART(predictors, param)
  τgrid <- grids$τgrid
  μgrid <- grids$μgrid
  dichotomous <- grids$dichotomous
  n <- grids$n
  p <- grids$p

  gamma0 <- initialize_gamma(outcomes, param)
  gammafit <- rep(gamma0, n)

  rh <- evaluate_pseudoresid(outcomes, gammafit)
  SMARTtrees <- SMARTboostTrees(param, gamma0, n, p, meanx, stdx)

  for (iter in 1:param$ntrees) {
    Gβ <- i <- μ <- τ <- β <- fi2 <- fit_one_tree(rh$r, rh$h, data$x, SMARTtrees$infeatures, μgrid, dichotomous, τgrid, param)

    tree <- list(i = i, μ = μ, τ = τ, β = β, fi2 = fi2)
    updateSMARTtrees(SMARTtrees, Gβ, tree, rh, iter)
    rh <- evaluate_pseudoresid(data$y, SMARTtrees$gammafit, param)
  }

  return(SMARTtrees)
}
