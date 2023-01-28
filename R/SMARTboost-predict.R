#' Predict from a `SMARTboost`
#'
#' @param object A `SMARTboost` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"numeric"` for numeric predictions.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @examples
#' train <- mtcars[1:20,]
#' test <- mtcars[21:32, -1]
#'
#' # Fit
#' mod <- SMARTboost(mpg ~ cyl + log(drat), train)
#'
#' # Predict, with preprocessing
#' predict(mod, test)
#'
#' @export
predict.SMARTboost <- function(object, new_data, type = "numeric", ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_SMARTboost_predict_types())
  predict_SMARTboost_bridge(type, object, forged$predictors)
}

valid_SMARTboost_predict_types <- function() {
  c("numeric")
}

# ------------------------------------------------------------------------------
# Bridge

predict_SMARTboost_bridge <- function(type, model, predictors) {
  predictors <- as.matrix(predictors)

  predict_function <- get_SMARTboost_predict_function(type)
  predictions <- predict_function(model, predictors)

  hardhat::validate_prediction_size(predictions, predictors)

  predictions
}

get_SMARTboost_predict_function <- function(type) {
  switch(
    type,
    numeric = predict_SMARTboost_numeric
  )
}

# ------------------------------------------------------------------------------
# Implementation

predict_SMARTboost_numeric <- function(model, predictors) {

  if(is.matrix(predictors)) {

    names(model$SMARTtrees$meanx)

    stopifnot("Column names does not match" = identical(colnames(predictors), names(model$SMARTtrees$meanx)))

    predictors <- predictors[, reorder(colnames(predictors), match(colnames(predictors), names(model$SMARTtrees$meanx)))]

    predictors <- (predictors - model$SMARTtrees$meanx)/model$SMARTtrees$stdx
    gammafit <- model$SMARTtrees$gamma0
    for(j in 1:length(model$SMARTtrees$trees)) {
      tree <- model$SMARTtrees$trees[[j]]
      gammafit <- gammafit + model$SMARTtrees$param$lambda*SMARTtreebuild(predictors, tree$i, tree$mu, tree$tau, tree$beta, model$SMARTtrees$param$sigmoid)
    }
  }else {
    stop("Input 'x' must be either a matrix or a vector.")
  }
  hardhat::spruce_numeric(gammafit[,1])
}
