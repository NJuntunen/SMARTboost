
updateG <- function(G0, g) {
  if(is.vector(G0)){
    p <- 1
    n <- length(G0)

    G <- matrix(nrow = n, ncol = p*2)
    for (i in 1:p) {
      G[, i] <- G0 * g
      G[, i+p] <- G0 * (1 - g)
    }
  }else{
    n <- nrow(G0)
    p <- ncol(G0)

    G <- matrix(nrow = n, ncol = p*2)
    for (i in 1:p) {
      G[, i] <- G0[,i] * g
      G[, i+p] <- G0[,i] * (1 - g)
    }
  }
  return(G)
}




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
  # dichotomous <- rep(FALSE, p)

  if (n > maxn) {
    ssi <- sample(1:n, maxn, replace = FALSE)
  } else {
    ssi <- 1:n
  }

  dt <- data.table::as.data.table(x)

  dichotomous <- purrr::map_vec(data.table::as.data.table(dt), function(col) {
    length(unique(col)) == 2
  })

  mgrid <- purrr::map2_df(dt, 1:p, function(col, i) {
    if (!dichotomous[i]) {
      mgrid_col <- quantile(col, (1/(npoints+1)):(1/(npoints+1)):(1-1/(npoints+1)), na.rm = TRUE)
      u <- unique(mgrid_col)
      lu <- length(u)
      if (lu <= 3) {
        mgrid_col[1:lu] <- u[1:lu]
        mgrid_col[(lu+1):npoints] <- quantile(unique(col), (1/(npoints+1-lu)):(1/(npoints+1-lu)):(1-1/(npoints+1-lu)), na.rm = TRUE)
      }
      return(mgrid_col)
    }
  })

  return (list(mgrid = as.matrix(mgrid), dichotomous = dichotomous))
}

loopfeatures <- function(r, h, G0, x, ifit, infeatures, mugrid, dichotomous, taugrid, param, var_epsilon,cl) {
  p <- ncol(x)
  outputarray <- tibble()
  ps <- 1:p

  if (param$subsampleshare_columns < 1) {
    psmall <- round(p*param$subsampleshare_columns)
    ps <- sample(ps, psmall)
  }

  t_fun <- function(i) {
    t <- list(r = r, h = h, G0 = G0, xi = x[, i], infeaturesfit = updateinfeatures(infeatures, i),
              dichotomous = dichotomous, mugridi = mugrid[, i], dichotomous_i = dichotomous[i], taugrid = taugrid, param = param, var_epsilon = var_epsilon)

    add_depth_cpp(t)
  }


  outputarray <- do.call(rbind, purrr::map(ps, t_fun))

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

Gfitbeta2 <- function(r, h, G0, xi, param, var_epsilon, infeaturesfit, dichotomous, muv, tau, dichotomous_i, G) {
  mu <- muv[1]
  tau <- max(tau, 0.2)  # Anything lower than 0.2 is still essentially linear, with very flat log-likelihood

  gL <- sigmoidf_cpp(xi, mu, tau, param$sigmoid, dichotomous = dichotomous_i)
  G <- updateG_allocated_cpp(G0, gL, G)

  loss <- fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit, dichotomous, mu, tau, dichotomous_i)

  return(loss$loss)
}

add_depth <- function(t) {

 lossmatrix <- matrix(Inf,length(t$taugrid),length(t$mugridi))
 n <- dim(t$G0)[1]
 p <- dim(t$G0)[2]
 G <- matrix(NA, n, 2*p)

  if(t$dichotomous_i==TRUE) {   # no optimization needed
    loss <- Gfitbeta_cpp(t$r,t$h,t$G0,t$xi,t$param,t$var_epsilon,t$infeaturesfit,t$dichotomous,c(0,0),t$dichotomous_i,G)
    tau <- 999.9
    mu <- 0
  } else {
    for (indexmu in 1:length(t$mugridi)) {
      for (indextau in 1:length(t$taugrid)) {
        lossmatrix[indextau,indexmu] <- Gfitbeta_cpp(t$r,t$h,t$G0,t$xi,t$param,t$var_epsilon,t$infeaturesfit,t$dichotomous,c(t$mugridi[indexmu],log(t$taugrid[indextau])),t$dichotomous_i,G)
        if (indextau > 1 && lossmatrix[indextau,indexmu] > lossmatrix[indextau-1,indexmu]) { break }   #  if loss increases, break loop over tau (reduces computation costs by some 25%)
      }
    }

    minindex <- which(lossmatrix == min(lossmatrix), arr.ind = TRUE)[1,]   # returns a Cartesian index
    loss <- lossmatrix[minindex[1], minindex[2]]
    tau <- t$taugrid[minindex[1]]
    mu <- t$mugridi[minindex[2]]

    # Optionally, further optimize over mu. Perhaps needed for highly nonlinear functions.
    if (t$param$optimizevs) {
      mu0 <- t$mu
      res <- nloptr(x0 = mu0,
                    eval_f = function(mu) Gfitbeta2(t$r, t$h, t$G0, t$xi, t$param, t$var_epsilon, t$infeaturesfit, t$dichotomous, mu, t$tau, t$dichotomous_i, G),
                    lb = rep(-Inf, length(mu0)),
                    ub = rep(Inf, length(mu0)),
                    opts = list("algorithm" = "NLOPT_LN_BOBYQA",
                                "xtol_rel" = t$param$xtolOptim,
                                "maxeval" = 100))
      loss <- res$objective
      mu <- res$solution[1]
    }
  }

  result <- tibble(loss = loss, tau = tau, mu = mu)
  return(result)
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

optimize_mutau <- function(r,h,G0,xi,param,var_epsilon,infeaturesfit,dichotomous,tau,dichotomous_i,mu0,T) {
  n <- dim(G0)[1]
  p <- dim(G0)[2]
  G <- matrix(NA,n,p*2)

  objective_function <- function(mu) {
    Gfitbeta2(r,h,G0,xi,param,var_epsilon,infeaturesfit,dichotomous,mu,tau,dichotomous_i,G)
  }

  options <- list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 100, xtol_rel = param$xtolOptim/(1+tau))
  optim_result <- nloptr::nloptr(x0 = mu0, eval_f = objective_function, opts = options)

  return(optim_result)
}


refineOptim <- function(r,h,G0,xi,infeaturesfit,dichotomous,mu0,dichotomous_i,tau0,param,var_epsilon,cl) {
  if (dichotomous_i) {
    gL  <- sigmoidf_cpp(xi,mu0,tau0,param$sigmoid,dichotomous=dichotomous_i)
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

    lossmatrix <- purrr::map(1:length(taugrid), function(index_tau) {
      res <- optimize_mutau(r,h,G0,xi,param,var_epsilon,infeaturesfit,dichotomous,taugrid[index_tau],dichotomous_i,mu0)
      c(res$objective, res$solution[1])
    })
    lossmatrix <- do.call(rbind, lossmatrix)

    minindex <- which.min(lossmatrix[,1])
    loss <- lossmatrix[minindex,1]
    tau <- taugrid[minindex]
    mu <- lossmatrix[minindex,2]
  }

  return(list(loss=loss,tau=tau,mu=mu))
}




fit_one_tree <- function(r, h, x, infeatures, mugrid, dichotomous, taugrid, param) {
  var_wr <- var(r)
  var_epsilon <- var_wr * (1 - param$R2p)

  n <- nrow(x)
  p <- ncol(x)
  G0 <- matrix(1, n, 1) # initialize G, the matrix of features
  loss0 <- Inf
  ifit <- vector()
  mufit <- vector()
  taufit <- vector()
  yfit0 <- rep(0, n)
  infeaturesfit <- infeatures
  fi2 <- rep(0, param$depth)

  subsamplesize <- round(n * param$subsamplesharevs)

  if (param$subsamplesharevs == 1) {
    ssi <- 1:n
  } else {
    ssi <- sample(1:n, subsamplesize, replace = FALSE) # sub-sample indexes. Sub-sample no reimmission
  }

  for (depth in 1:param$depth) { #  NB must extend G for this to be viable
    # variable selection

    if (param$subsamplesharevs == 1) {
      for(i in 1:100){
        outputarray <- loopfeatures_cpp(r, h, G0, x, ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, var_epsilon) # loops over all variables
      }
    } else { # Variable selection using a random sub-set of the sample. All the sample is then used in refinement.
      if (length(h) == 1) {
        outputarray <- loopfeatures(r[ssi], h[ssi], G0[ssi,], x[ssi,], ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, var_epsilon) # loops over all variables
      } else {
        outputarray <- loopfeatures(r[ssi], h, G0[ssi,], x[ssi,], ifit, infeaturesfit, mugrid, dichotomous, taugrid, param, var_epsilon) # loops over all variables
      }
    }
    #outputarray <- do.call(rbind, outputarray)
    i <- which(outputarray[, 1] == min(outputarray[, 1]), arr.ind = TRUE)[1] # outputarray[,1] is loss (minus log marginal likelihood) vector
    tau0 <- outputarray[i, 2]
    mu0 <- outputarray[i, 3]

    infeaturesfit <- updateinfeatures(infeaturesfit, i)

    # refine optimization, after variable selection
    if (param$subsamplesharevs < 1 && param$subsamplefinalbeta == TRUE) {
      if (length(h) == 1) {
        loss <- refineOptim(r[ssi], h[ssi], G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, var_epsilon,cl)
      } else {
        loss <- refineOptim(r[ssi], h, G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, var_epsilon,cl)
      }
    } else {
      loss <- refineOptim(r, h, G0, x[ssi, i], infeaturesfit, dichotomous, mu0, dichotomous[i], tau0, param, var_epsilon,cl)
    }
    tau <- loss$tau
    mu <- loss$mu
    loss <- loss$loss

    gL <- sigmoidf_cpp(x[, i], mu, tau, param$sigmoid, dichotomous = dichotomous[i])
    G <- updateG_allocated_cpp(G0, gL, matrix(NA, n, 2^depth))

    loss <- fitbeta_cpp(r, G, var_epsilon, param, infeaturesfit, dichotomous, mu, tau, dichotomous[i])
    yfit <- loss$Gbeta
    beta <- loss$beta
    loss <- loss$loss

    fi2[depth] <- (sum(yfit^2) - sum(yfit0^2)) / n

    G0 <- G
    loss0 <- loss
    yfit0 <- yfit
    ifit <- c(ifit, i)
    mufit <- c(mufit, mu)
    taufit <- c(taufit, tau)
    betafit <- beta

  }


  return(list(yfit0=yfit0,ifit=ifit,mufit=mufit,taufit=taufit,betafit=betafit,fi2=fi2))

}

SMARTtreebuild <- function(x, ij, muj, tauj, betaj, sigmoid) {
  n <- nrow(x)
  p <- ncol(x)
  depth <- length(ij)
  G <- rep(1, n)

  for (d in 1:depth) {
    i <- ij[d]
    mu <- muj[d]
    tau <- tauj[d]
    gL <- sigmoidf_cpp(x[,i], mu, tau, sigmoid)
    G <- updateG(G, gL)
  }

  gammafit <- G * betaj[,1]
  return(gammafit)
}














