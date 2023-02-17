fitbeta <- function(r, h, G, param, var_epsilon, infeaturesfit, dichotomous, mu, tau, dichotomous_i) {

  var_r <- var_epsilon/(1-param$R2p)  # var_epsilon is computed on the pseudo-residuals
  n <- nrow(G)
  p <- ncol(G)

  GGh <- t(G) %*% G

  for (i in 1:p) {
    GGh[i,i] <- max(GGh[i,i], max(diag(GGh)*0.00001))
  }

  Pb <- (sum(diag(GGh)) / (n*var_r*param$R2p))
  beta <- matrix(0,nrow=p,ncol=1)

  tryCatch({beta <- solve(GGh + var_r*param$loglikdivide*Pb, t(G) %*% r)}, warning = function(w) {
    while(class(w)=="SingularException") {
      Pb <- Pb*2.01
      beta <- solve(GGh + var_r*param$loglikdivide*Pb, t(G) %*% r)
    }
  })

  Gbeta <- G %*% beta

  loglik <- -0.5*( (sum((r - Gbeta)^2)/var_r) )/param$loglikdivide

  logpdfbeta <- -0.5*( p*log(2*pi) - p*log(Pb) + Pb*(t(beta) %*% beta) )

  if (dichotomous_i) {
    logpdfmu <- 0
    logpdftau <- 0
  } else if (param$sharptree==TRUE) {
    logpdfmu <- lnpmu(mu,param$varmu,param$dofmu)
    logpdftau <- 0
  } else {
    logpdfmu <- lnpmu(mu,param$varmu,param$dofmu)
    logpdftau <- lnptau(tau,param$meanlntau,param$varlntau,param$doflntau,param$depth)
  }

  loss <- -( loglik + logpdfbeta + logpdftau + logpdfmu)

  return(list(loss=loss, Gbeta=Gbeta, beta = beta))
}


Gfitbeta <- function(r, h, G0, xi, param, var_epsilon, infeaturesfit, dichotomous, muLogTau, dichotomous_i, G) {

  mu <- muLogTau[1]
  tau <- exp(muLogTau[2])
  tau <- max(tau, 0.2) # Anything lower than 0.2 is still essentially linear, with very flat log-likelihood

  gL <- sigmoidf(xi, mu, tau, param$sigmoid, dichotomous = dichotomous_i)
  G <- updateG_allocated(G0, gL, G)

  result <- fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit, dichotomous, mu, tau, dichotomous_i)

  return(result$loss)
}
