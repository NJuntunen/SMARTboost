#' Wrapper to add `SMARTboost` engine to the parsnip `boost_tree` model
#' specification
#'
#' @return NULL
#' @export
add_parsnip_SMARTboost <- function() {
  parsnip::set_new_model("SMARTboost")
  parsnip::set_model_mode(model = "SMARTboost", mode = "regression")

  parsnip::set_model_engine(
    "SMARTboost",
    mode = "regression",
    eng = "SMARTboost"
  )

  parsnip::set_dependency("SMARTboost", eng = "SMARTboost", pkg = "SMARTboost")

  parsnip::set_fit(
    model = "SMARTboost",
    eng = "SMARTboost",
    mode = "regression",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "SMARTboost", fun = "SMARTboost_fit"),
      defaults = list()
    )
  )

  parsnip::set_encoding(
    model = "SMARTboost",
    eng = "SMARTboost",
    mode = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  make_class_info <- function(type) {
    list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args =
        list(
          object = quote(object$fit),
          new_data = quote(new_data),
          type = type
        )
    )
  }

  parsnip::set_model_arg(
    model = "SMARTboost",
    eng = "SMARTboost",
    parsnip = "tree_depth",
    original = "depth",
    func = list(pkg = "dials", fun = "tree_depth"),
    has_submodel = TRUE
  )

  parsnip::set_model_arg(
    model = "SMARTboost",
    eng = "SMARTboost",
    parsnip = "trees",
    original = "ntrees",
    func = list(pkg = "dials", fun = "trees"),
    has_submodel = FALSE
  )

  parsnip::set_model_arg(
    model = "SMARTboost",
    eng = "SMARTboost",
    parsnip = "learn_rate",
    original = "lambda",
    func = list(pkg = "dials", fun = "learn_rate"),
    has_submodel = FALSE
  )

  parsnip::set_model_arg(
    model = "SMARTboost",
    eng = "SMARTboost",
    parsnip = "mtry",
    original = "subsampleshare_columns",
    func = list(pkg = "dials", fun = "mtry"),
    has_submodel = FALSE
  )
}

#' @export
SMARTboost <- function(mode = "regression", tree_depth = NULL, trees = NULL, learn_rate = NULL,
                       mtry = NULL) {

  if (!requireNamespace("parsnip", quietly = TRUE))
    rlang::abort("Package \"parsnip\" needed for this function to work. Please install it.")

  # Capture the arguments in quosures
  args <- list(
    tree_depth = rlang::enquo(tree_depth),
    trees = rlang::enquo(trees),
    learn_rate = rlang::enquo(learn_rate),
    mtry = rlang::enquo(mtry)
  )

  # Save some empty slots for future parts of the specification
  out <- list(args = args, eng_args = NULL,
              mode = mode, method = NULL, engine = NULL)

  # set classes in the correct order
  class(out) <- parsnip::make_classes("SMARTboost")
  out
}

#' @export
#' @importFrom stats update
update.SMARTboost <- function(object, parameters = NULL, depth = NULL, ntrees = NULL, lambda = NULL,
                              subsampleshare_columns = NULL, ...) {
  rlang::check_installed("parsnip")
  eng_args <- parsnip::update_engine_parameters(object$eng_args, fresh=TRUE, ...)
  args <- list(
    tree_depth = rlang::enquo(tree_depth),
    trees = rlang::enquo(trees),
    learn_rate = rlang::enquo(learn_rate),
    mtry = rlang::enquo(mtry)
  )
  args <- parsnip::update_main_parameters(args, parameters)
  parsnip::new_model_spec(
    "SMARTboost",
    args = args,
    eng_args = eng_args,
    mode = object$mode,
    method = NULL,
    engine = object$engine
  )
}



.onLoad <- function(libname, pkgname) {

  add_parsnip_SMARTboost()
}

