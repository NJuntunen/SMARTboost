SMARTbst <- function(data, param) {
  # initialize SMARTtrees
  data <- preparedataSMART(data, param)
  meanx <- data$meanx
  stdx <- data$stdx
  τgrid <- data$τgrid
  μgrid <- data$μgrid
  dichotomous <- data$dichotomous
  n <- data$n
  p <- data$p
  
  gamma0 <- initialize_gamma(data, param)
  gammafit <- rep(gamma0, n)
  
  rh <- evaluate_pseudoresid(data$y, gammafit, param)
  SMARTtrees <- SMARTboostTrees(param, gamma0, n, p, meanx, stdx)
  
  for (iter in 1:param$ntrees) {
    Gβ <- i <- μ <- τ <- β <- fi2 <- fit_one_tree(rh$r, rh$h, data$x, SMARTtrees$infeatures, μgrid, dichotomous, τgrid, param)
    
    tree <- list(i = i, μ = μ, τ = τ, β = β, fi2 = fi2)
    updateSMARTtrees(SMARTtrees, Gβ, tree, rh, iter)
    rh <- evaluate_pseudoresid(data$y, SMARTtrees$gammafit, param)
  }
  
  return(SMARTtrees)
}
