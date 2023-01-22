new_SMARTboost <- function(SMARTtrees, blueprint) {
  hardhat::new_model(SMARTtrees = SMARTtrees, blueprint = blueprint, class = "SMARTboost")
}
