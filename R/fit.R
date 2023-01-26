gridvectortau <- function(meanlntau, taugridpoints, sharptree = FALSE) {
  stopifnot("taugridpoints must be between 3 and 5" = 1 <= taugridpoints & taugridpoints <= 5)
  if (sharptree) {
    taugrid <- c(100.0)
  } else {
    if (taugridpoints == 1) {
      taugrid <- c(5.0)
    } else if (taugridpoints == 2) {
      taugrid <- c(1.0, 9.0)
    } else if (taugridpoints == 3) {
      taugrid <- c(1.0, 3.0, 9.0)
    } else if (taugridpoints == 4) {
      taugrid <- c(1.0, 3.0, 9.0, 27.0)
    } else if (taugridpoints == 5) {
      taugrid <- c(1.0, 3.0, 9.0, 27.0, 50.0)
    }
  }
  return(taugrid)
}



gridmatrixmu <- function(x, npoints, tol = 0.005, maxiter = 100, fuzzy = FALSE, maxn = 100000) {
  n <- nrow(x)
  p <- ncol(x)
  stopifnot("npoints cannot be larger than n" = npoints < n)
  mgrid <- matrix(NA, nrow = npoints, ncol = p)
  dichotomous <- rep(FALSE, p)

  if (n > maxn) {
    ssi <- sample(1:n, maxn, replace = FALSE)
  } else {
    ssi <- 1:n
  }

  for (i in 1:p) {
    dichotomous[i] <- length(unique(x[,i])) == 2
    if (!dichotomous[i]) {
      mgrid[,i] <- quantile(x[,i], (1/(npoints+1)):(1/(npoints+1)):(1-1/(npoints+1)))
      u <- unique(mgrid[,i])
      lu <- length(u)
      if (lu <= 3) {
        mgrid[1:lu,i] <- u[1:lu]
        mgrid[(lu+1):npoints,i] <- quantile(unique(x[,i]), (1/(npoints+1-lu)):(1/(npoints+1-lu)):(1-1/(npoints+1-lu)))
      }
    }
  }

  return (list(mgrid = mgrid, dichotomous = dichotomous))
}

loopfeatures <- function(r, h, G0, x, ifit, infeatures, mugrid, dichotomous, taugrid, param, vareps) {
  p <- ncol(x)
  outputarray <- matrix(rep(Inf, p*3), nrow = p, ncol = 3)
  ps <- 1:p

  if (param$subsampleshare_columns < 1) {
    psmall <- round(p*param$subsampleshare_columns)
    ps <- sample(ps, psmall)
  }

  for (i in ps) {
    t <- list(r = r, h = h, G0 = G0, xi = x[,i], infeaturesfit = updateinfeatures(infeatures, i),
              dichotomous = dichotomous, mugridi = mugrid[,i], dichotomous_i = dichotomous[i], taugrid = taugrid, param = param, vareps = vareps)
    outputarray[i, ] <- add_depth(t)
  }
  return(outputarray)
}

sigmoidf <- function(x, mu, tau, sigmoid, dichotomous = FALSE) {
  if(dichotomous){   # return 0 if x<=0 and x>1 otherwise. x is assumed de-meaned
    g <- ifelse(x > 0, 1, 0)
    return(g)
  }

  if(sigmoid == "sigmoidsqrt"){
    g <- 0.5 + 0.5 * ( 0.5 * tau * (x - mu) / sqrt( (1 + (0.5 * tau * (x - mu))^2) ) )
  }else if(sigmoid == "sigmoidlogistic"){
    g <- 1 - 1/(1 + exp(tau * (x - mu)))
  }

  return(g)
}

updateG_allocated <- function(G0, g, G) {
  p <- ncol(G0)
  for (i in 1:p) {
    G[,i] <- G0[,i] * g
    G[,i+p] <- G0[,i] * (1 - g)   #  =  G0[:,i] - G[:,i] is not faster
  }
  return(G)
}

lnpmu <- function(mu, varmu, dofmu) {
  s <- sqrt(varmu)
  lnp <- sum(dnorm(mu, mean = 0, sd = s, log = TRUE))
  return(lnp)
}

logpdft <- function(x, m, s, v) {
  z <- (x - m)/s
  logpdfz <- -0.5723649429247001 + lgamma((v + 1)/2) - lgamma(v/2) - 0.5*log(v) - 0.5*(1 + v)*log(1 + (z^2)/v)
  return(logpdfz - log(s))
}


lnptau <- function(tau, meanlntau, varlntau, doflntau, depth) {
  s <- sqrt(varlntau/depth)
  lnp <- sum(logpdft(log(tau), meanlntau, s, doflntau))
  return(lnp)
}


fitbeta <- function(r, h, G, param, var_epsilon, infeaturesfit, dichotomous, mu, tau, dichotomous_i) {

  var_r <- var_epsilon/(1-param$R2p)  # var_epsilon is computed on the pseudo-residuals
  n <- nrow(G)
  p <- ncol(G)

  GGh <- t(G) %*% G

  Pb <- (sum(diag(GGh))/(n*var_r*param$R2p))*diag(p)

  beta <- matrix(0,nrow=p,ncol=1)

  tryCatch({beta <- solve(GGh + var_epsilon*param$loglikdivide*Pb, t(G) %*% r)}, warning = function(w) {
    while(class(w)=="SingularException") {
      Pb <- Pb*2.01
      beta <- solve(GGh + var_epsilon*param$loglikdivide*Pb, t(G) %*% r)
    }
  })

  Gbeta <- G %*% beta

  loglik <- -0.5*( (sum((r - Gbeta)^2)/var_epsilon) )/param$loglikdivide

  logpdfbeta <- -0.5*( p*log(2*pi) - p*log(Pb[1,1]) + Pb[1,1]*(t(beta) %*% beta) )

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

Gfitbeta <- function(r, h, G0, xi, param, varEpsilon, infeaturesfit, dichotomous, muLogTau, dichotomous_i, G) {

  mu <- muLogTau[1]
  tau <- exp(muLogTau[2])
  tau <- max(tau, 0.2) # Anything lower than 0.2 is still essentially linear, with very flat log-likelihood

  gL <- sigmoidf(xi, mu, tau, param$sigmoid, dichotomous = dichotomous_i)
  G <- updateG_allocated(G0, gL, G)

  result <- fitbeta(r, h, G, param, varEpsilon, infeaturesfit, dichotomous, mu, tau, dichotomous_i)

  return(result$loss)
}

Gfitbeta2 <- function(r, h, G0, xi, param, vareps, infeaturesfit, dichotomous, muv, tau, dichotomous_i, G) {
  mu <- muv[1]
  tau <- max(tau, 0.2)  # Anything lower than 0.2 is still essentially linear, with very flat log-likelihood

  gL <- sigmoidf(xi, mu, tau, param$sigmoid, dichotomous = dichotomous_i)
  G <- updateG_allocated(G0, gL, G)

  loss <- fitbeta(r, h, G, param, vareps, infeaturesfit, dichotomous, mu, tau, dichotomous_i)

  return(loss)
}

#TODO : tama funktio ei toimi oikein, kysy chatGPTltä:
#Write this julia optimization code using R and nloptr package:  if t.param.optimizevs==true
#μ0   = [μ]
#res  = Optim.optimize( μ -> Gfitβ2(t.r,t.h,t.G0,t.xi,t.param,t.varϵ,t.infeaturesfit,t.dichotomous,μ,τ,t.dichotomous_i,G),μ0,Optim.BFGS(linesearch = LineSearches.BackTracking()), Optim.Options(iterations = 100,x_tol = t.param.xtolOptim))
#loss = res.minimum
#μ    = res.minimizer[1]
#end
#add_depth <- function(t) {
#  T <- typeof(t$varϵ)
#  lossmatrix <- matrix(Inf,length(t$τgrid),length(t$μgridi))

#  n <- dim(t$G0)[1]
#  p <- dim(t$G0)[2]
#  G   <- matrix(numeric(n*(2*p)), n, 2*p)

  if(t$dichotomous_i==TRUE) {   # no optimization needed
    loss <- Gfitβ(t$r,t$h,t$G0,t$xi,t$param,t$varϵ,t$infeaturesfit,t$dichotomous,c(0,0),t$dichotomous_i,G)
    τ <- 999.9
    μ <- 0
  } else {
    for (indexμ in 1:length(t$μgridi)) {
      for (indexτ in 1:length(t$τgrid)) {
        lossmatrix[indexτ,indexμ] <- Gfitβ(t$r,t$h,t$G0,t$xi,t$param,t$varϵ,t$infeaturesfit,t$dichotomous,c(t$μgridi[indexμ],log(t$τgrid[indexτ])),t$dichotomous_i,G)
        if (indexτ > 1 && lossmatrix[indexτ,indexμ] > lossmatrix[indexτ-1,indexμ]) { break }   #  if loss increases, break loop over tau (reduces computation costs by some 25%)
      }
    }

    minindex <- which(lossmatrix == min(lossmatrix), arr.ind = TRUE)   # returns a Cartesian index
    loss <- lossmatrix[minindex]
    tau <- t$taugrid[minindex[1]]
    mu <- t$mugridi[minindex[2]]

    # Optionally, further optimize over mu. Perhaps needed for highly nonlinear functions.
    if (t$param$optimizevs == TRUE) {
      mu0 <- mu
      res <- optim(
        function(mu) Gfitbeta2(t$r, t$h, t$G0, t$xi, t$param, t$vareps, t$infeaturesfit, t$dichotomous, mu, tau, t$dichotomous_i, G),
        par = mu0,
        method = "BFGS",
        control = list(iter.max = 100, reltol = t$param$xtolOptim)
      )
      loss <- res$value
      mu <- res$par
    }
  }
  return(c(loss, tau, mu))
}

updateinfeatures <- function(infeatures, ifit) {
  x <- infeatures
  for (i in ifit) {
    x[i] <- TRUE
  }
  return (x)
}

preparedataSMART <- function(y, data) {

  meanx <- apply(data, 2, mean)
  stdxL2 <- apply(data, 2, sd)
  stdx <- stdxL2

  for (i in 1:ncol(data)) {
    if (length(unique(data[,i])) > 2) { # with 0-1 data, stdx = 0 when computed on non-zero data.
      meanx[i] <- median(data[,i]) # median rather than mean
      xi <- data[,i]
      s_median <- 1.42 * median(abs(xi - meanx[i])) # use of robust measure of std, not std, and computed only on non-zero values
      stdx[i] <- ifelse(s_median == 0, stdxL2[i], min(s_median, stdxL2[i]))
    }
  }

  data_standardized <- list(y = y, x = (data - meanx) / stdx)

  return(list(data_standardized = data_standardized, meanx = meanx, stdx = stdx))
}

preparegridsSMART <- function(data, param) {
  taugrid <- gridvectortau(param$meanlntau, param$taugridpoints, sharptree = param$sharptree)
  gridmatrixmu <- gridmatrixmu(data, param$mugridpoints)
  mugrid <- gridmatrixmu[[1]]
  dichotomous <- gridmatrixmu[[2]]
  n <- nrow(data)
  p <- ncol(data)

  return(list(taugrid = taugrid, mugrid = mugrid, dichotomous = dichotomous, n = n, p = p))
}

optimize_mutau <- function(r,h,G0,xi,param,vareps,infeaturesfit,dichotomous,tau,dichotomous_i,mu0,T) {
  n <- dim(G0)[1]
  p <- dim(G0)[2]
  G <- matrix(NA,n,p*2)

  objective_function <- function(mu) {
    Gfitbeta2(r,h,G0,xi,param,vareps,infeaturesfit,dichotomous,mu,tau,dichotomous_i,G)
  }

  options <- list(algorithm = "NLOPT_LD_LBFGS", maxeval = 100, xtol_rel = param$xtolOptim/(1+tau>10))
  optim_result <- nloptr::nloptr(x0 = mu0, eval_f = objective_function, opts = options)

  return(optim_result)
}


refineOptim <- function(r,h,G0,xi,infeaturesfit,dichotomous,mu0,dichotomous_i,tau0,param,var_epsilon) {
  if (dichotomous_i) {
    gL  <- sigmoidf(xi,mu0,tau0,param$sigmoid,dichotomous=dichotomous_i)
    n <- dim(G0)[1]
    p <- dim(G0)[2]
    G   <- matrix(NA, n, p*2)
    G <- updateG_allocated(G0,gL,G)
    loss <- Inf
    tau <- tau0
    mu <- mu0
  } else {
    if (param$sharptree == TRUE) {
      taugrid <- tau0
    } else if (param$taugridpoints == 1) {
      taugrid <- exp(log(tau0) + seq(-2.7, 2.7, by = 0.3))
    } else if (param$taugridpoints == 2) {
      taugrid <- exp(log(tau0) + seq(-1.8, 1.8, by = 0.3))
    } else {
      if (tau0 < 8.0) {
        taugrid <- exp(log(tau0) + seq(-0.9, 0.9, by = 0.3))
      } else {
        taugrid <- exp(log(tau0) + seq(-0.9, 1.8, by = 0.3))
      }
    }

    lossmatrix <- matrix(NA,length(taugrid),2)
    for (index_tau in 1:length(taugrid)) {
      res <- optimize_mutau(r,h,G0,xi,param,var_epsilon,infeaturesfit,dichotomous,taugrid[index_tau],dichotomous_i,mu0)
      lossmatrix[index_tau,1] <- res$objective
      lossmatrix[index_tau,2] <- res$solution[1]
    }

    minindex <- which.min(lossmatrix[,1])
    loss <- lossmatrix[minindex,1]
    tau <- taugrid[minindex]
    mu <- lossmatrix[minindex,2]
  }

  return(list(loss=loss,tau=tau,mu=mu))
}




fit_one_tree <- function(r, h, x, infeatures, mugrid, dichotomous, taugrid, param) {
  var_wr <- var(r)
  vareps <- var_wr * (1 - param$R2p)

  n <- nrow(x)
  p <- ncol(x)
  G0 <- matrix(1, n, 1) # initialize G, the matrix of features
  loss0 <- Inf

  yfit0 <- matrix(0, n, 1)
  ifit <- integer()
  mufit <- numeric()
  taufit <- numeric()
  infeaturesfit <- infeatures
  fi2 <- matrix(0, param$depth, 1)
  betafit <- numeric()

  subsamplesize <- round(n * param$subsamplesharevs)

  if (param$subsamplesharevs == 1) {
    ssi <- 1:n
  } else {
    ssi <- sample(1:n, subsamplesize, replace = FALSE) # sub-sample indexes. Sub-sample no reimmission
  }

  for (depth in 1:param$depth) { #  NB must extend G for this to be viable
    # variable selection
    if (param$subsamplesharevs == 1) {
      outputarray <- loopfeatures(r, h, G0, x, ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, vareps) # loops over all variables
    } else { # Variable selection using a random sub-set of the sample. All the sample is then used in refinement.
      if (length(h) == 1) {
        outputarray <- loopfeatures(r[ssi], h[ssi], G0[ssi,], x[ssi,], ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, vareps) # loops over all variables
      } else {
        outputarray <- loopfeatures(r[ssi], h, G0[ssi,], x[ssi,], ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, vareps) # loops over all variables
      }
    }

    i <- which.min(outputarray[, 1]) # outputarray[,1] is loss (minus log marginal likelihood) vector
    tau0 <- outputarray[i, 2]
    mu0 <- outputarray[i, 3]
    print(outputarray)

    infeaturesfit <- updateinfeatures(infeaturesfit, i)

    # refine optimization, after variable selection
    if (param$subsamplesharevs < 1 && param$subsamplefinalbeta == TRUE) {
      if (length(h) == 1) {
        loss <- refineOptim(r[ssi], h[ssi], G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, vareps)
      } else {
        loss <- refineOptim(r[ssi], h, G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, vareps)
      }
    } else {
      loss <- refineOptim(r, h, G0, x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, vareps)
    }

    gL <- sigmoidf(x[, i], mu, tau, param$sigmoid, dichotomous = dichotomous[i])
    G <- updateG_allocated(G0, gL, matrix(NA, n, 2^depth))

    loss <- fitbeta(r, h, G, param, vareps, infeaturesfit, dichotomous, mu, tau, dichotomous[i])
    yfit <- loss$Gbeta
    beta <- loss$beta
    loss <- loss$loss

    fi2[depth] <- (sum(yfit^2) - sum(yfit0^2)) / n

    G0 <- G
    loss0 <- loss
    yfit0 <- yfit
    betafit <- beta
    ifit <- c(ifit, i)
    mufit <- c(mufit, mu)
    taufit <- c(taufit, tau)
    betafit <- beta

  }

  return(list(yfit0=yfit0,ifit=ifit,mufit=mufit,taufit=taufit,betafit=betafit,fi2=fi2))

}

