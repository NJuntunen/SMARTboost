#' Fit a `SMARTboost`
#'
#' `SMARTboost()` fits a model.
#'
#' @param x Depending on the context:
#'
#'   * A __data frame__ of predictors.
#'   * A __matrix__ of predictors.
#'   * A __recipe__ specifying a set of preprocessing steps
#'     created from [recipes::recipe()].
#'
#' @param y When `x` is a __data frame__ or __matrix__, `y` is the outcome
#' specified as:
#'
#'   * A __data frame__ with 1 numeric column.
#'   * A __matrix__ with 1 numeric column.
#'   * A numeric __vector__.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both the predictors and the outcome.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @return
#'
#' A `SMARTboost` object.
#'
#' @examples
#' predictors <- mtcars[, -1]
#' outcome <- mtcars[, 1]
#'
#' # XY interface
#' mod <- SMARTboost(predictors, outcome)
#'
#' # Formula interface
#' mod2 <- SMARTboost(mpg ~ ., mtcars)
#'
#' # Recipes interface
#' library(recipes)
#' rec <- recipe(mpg ~ ., mtcars)
#' rec <- step_log(rec, disp)
#' mod3 <- SMARTboost(rec, mtcars)
#'
#' @export
SMARTboost <- function(x, depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  UseMethod("SMARTboost")
}

#' @export
#' @rdname SMARTboost
SMARTboost.default <- function(x, depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  stop("`SMARTboost()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname SMARTboost
SMARTboost.data.frame <- function(x, y,depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  processed <- hardhat::mold(x, y)
  SMARTboost_bridge(processed, depth, ntrees, lambda, subsampleshare_columns, ...)
}

# XY method - matrix

#' @export
#' @rdname SMARTboost
SMARTboost.matrix <- function(x, y,depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  processed <- hardhat::mold(x, y)
  SMARTboost_bridge(processed, depth, ntrees, lambda, subsampleshare_columns, ...)
}

# Formula method

#' @export
#' @rdname SMARTboost
SMARTboost.formula <- function(formula, data,depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  processed <- hardhat::mold(formula, data)
  SMARTboost_bridge(processed, depth, ntrees, lambda, subsampleshare_columns, ...)
}

# Recipe method

#' @export
#' @rdname SMARTboost
SMARTboost.recipe <- function(x, data,depth = 4, ntrees = 100, lambda = 0.2, subsampleshare_columns = 1, ...) {
  processed <- hardhat::mold(x, data)
  SMARTboost_bridge(processed, depth, ntrees, lambda, subsampleshare_columns, ...)
}

# ------------------------------------------------------------------------------
# Bridge

SMARTboost_bridge <- function(processed, depth = NULL, ntrees = NULL, lambda = NULL, subsampleshare_columns = NULL, ...) {

  param <- SMARTparam(depth = depth,
                      ntrees = ntrees,
                      lambda = lambda,
                      subsampleshare_columns = subsampleshare_columns,
                      ... = ...)

  predictors <- processed$predictors %>%
    as.matrix()
  outcomes <- processed$outcomes[[1]]

  date <- processed$extras$roles$date

  # initialize SMARTtrees
  data <- preparedataSMART(outcomes, predictors)
  meanx <- data$meanx
  stdx <- data$stdx
  data <- data$data_standardized

  grids <- preparegridsSMART(data$x, param)
  taugrid <- grids$taugrid
  mugrid <- grids$mugrid
  dichotomous <- grids$dichotomous
  n <- grids$n
  p <- grids$p

  gamma0 <- initialize_gamma(data$y, param)
  gammafit <- rep(gamma0, n)

  rh <- evaluate_pseudoresid(data$y, gammafit)
  SMARTtrees <- SMARTboostTrees(param, gamma0, n, p, meanx, stdx)

  for (iter in 1:param$ntrees) {
    out <- fit_one_tree(rh$r, rh$h, predictors,SMARTtrees$infeatures,mugrid, dichotomous, taugrid, param)

    tree <- list(i = out$ifit, mu = out$mufit, tau = out$taufit, beta = out$betafit, fi2 = out$fi2)
    SMARTtrees <- updateSMARTtrees(SMARTtrees, out$yfit0, tree, rh, iter)
    rh <- evaluate_pseudoresid(data$y, SMARTtrees$gammafit)
  }

  new_SMARTboost(
    SMARTtrees = SMARTtrees,
    blueprint = processed$blueprint
  )

}
