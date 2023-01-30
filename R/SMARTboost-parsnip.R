#' Wrapper to add `SMARTboost` engine to the parsnip `boost_tree` model
#' specification
#'
#' @return NULL
#' @export
add_boost_tree_SMARTboost <- function() {
  parsnip::set_model_engine("boost_tree", mode = "regression", eng = "SMARTboost")
  parsnip::set_dependency("boost_tree", eng = "SMARTboost", pkg = "SMARTboost")

  parsnip::set_fit(
    model = "boost_tree",
    eng = "SMARTboost",
    mode = "regression",
    value = list(
      interface = "data.frame",
      protect = c("x", "y"),
      func = c(pkg = "SMARTboost", fun = "SMARTboost"),
      defaults = list()
    )
  )

  parsnip::set_encoding(
    model = "boost_tree",
    mode = "regression",
    eng = "SMARTboost",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  parsnip::set_pred(
    model = "boost_tree",
    eng = "SMARTboost",
    mode = "regression",
    type = "numeric",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(pkg = "SMARTboost", fun = "predict.SMARTboost"),
      args = list(
        object = quote(object),
        new_data = quote(new_data)
      )
    )
  )

  # model args ----------------------------------------------------
  parsnip::set_model_arg(
    model = "boost_tree",
    eng = "SMARTboost",
    parsnip = "tree_depth",
    original = "depth",
    func = list(pkg = "dials", fun = "tree_depth"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "boost_tree",
    eng = "SMARTboost",
    parsnip = "trees",
    original = "ntrees",
    func = list(pkg = "dials", fun = "trees"),
    has_submodel = TRUE
  )
  parsnip::set_model_arg(
    model = "boost_tree",
    eng = "SMARTboost",
    parsnip = "learn_rate",
    original = "lambda",
    func = list(pkg = "dials", fun = "learn_rate"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "boost_tree",
    eng = "SMARTboost",
    parsnip = "mtry",
    original = "subsampleshare_columns",
    func = list(pkg = "dials", fun = "mtry"),
    has_submodel = FALSE
  )
}
