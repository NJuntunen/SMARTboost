gridvectorτ <- function(meanlnτ, τgridpoints, sharptree = FALSE) {
  stopifnot("τgridpoints must be between 3 and 5" = 1 <= τgridpoints & τgridpoints <= 5)
  if (sharptree) {
    τgrid <- c(100.0)
  } else {
    if (τgridpoints == 1) {
      τgrid <- c(5.0)
    } else if (τgridpoints == 2) {
      τgrid <- c(1.0, 9.0)
    } else if (τgridpoints == 3) {
      τgrid <- c(1.0, 3.0, 9.0)
    } else if (τgridpoints == 4) {
      τgrid <- c(1.0, 3.0, 9.0, 27.0)
    } else if (τgridpoints == 5) {
      τgrid <- c(1.0, 3.0, 9.0, 27.0, 50.0)
    }
  }
  return(τgrid)
}



gridmatrixμ <- function(x, npoints, tol = 0.005, maxiter = 100, fuzzy = FALSE, maxn = 100000) {
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

loopfeatures <- function(r, h, G0, x, ifit, infeatures, μgrid, dichotomous, τgrid, param, varϵ) {
  p <- ncol(x)
  outputarray <- matrix(Inf, p, 3)

  ps <- 1:p

  if (param$subsampleshare_columns < 1) {
    psmall <- round(p * param$subsampleshare_columns)
    ps <- sample(ps, psmall)
  }

  for (i in ps) {
    t <- list(r = r, h = h, G0 = G0, xi = x[, i], infeaturesfit = updateinfeatures(infeatures, i), dichotomous = dichotomous, μgridi = μgrid[, i], dichotomous_i = dichotomous[i], τgrid = τgrid, param = param, varϵ = varϵ)
    outputarray[i, ] <- add_depth(t)
  }

  outputarray
}

# TODO check that correct
add_depth <- function(t) {
  lossmatrix <- matrix(Inf, length(t$τgrid), length(t$μgridi))

  n <- nrow(t$G0)
  p <- ncol(t$G0)
  G <- matrix(NA, n, 2*p)

  if (t$dichotomous_i == TRUE) {
    # no optimization needed
    loss <- Gfitβ(t$r, t$h, t$G0, t$xi, t$param, t$varϵ, t$infeaturesfit, t$dichotomous, c(0, 0), t$dichotomous_i, G)
    τ <- 999.9
    μ <- 0
  } else {
    for (indexμ in 1:length(t$μgridi)) {
      for (indexτ in 1:length(t$τgrid)) {
        lossmatrix[indexτ, indexμ] <- Gfitβ(t$r, t$h, t$G0, t$xi, t$param, t$varϵ, t$infeaturesfit, t$dichotomous, c(t$μgridi[indexμ], log(t$τgrid[indexτ])), t$dichotomous_i, G)
        if (indexτ > 1 && lossmatrix[indexτ, indexμ] > lossmatrix[indexτ - 1, indexμ]) {
          break
        }
      }
    }

    minindex <- which.min(lossmatrix)
    loss <- lossmatrix[minindex]
    τ <- t$τgrid[minindex[1]]
    μ <- t$μgridi[minindex[2]]

    # Optionally, further optimize over μ. Perhaps needed for highly nonlinear functions.
    if (t$param$optimizevs == TRUE) {
      μ0 <- μ
      res <- optim(
        function(mu) Gfitβ2(t$r, t$h, t$G0, t$xi, t$param, t$varϵ, t$infeaturesfit, t$dichotomous, mu, τ, t$dichotomous_i, G),
        par = μ0,
        method = "BFGS",
        control = list(iter.max = 100, reltol = t$param$xtolOptim)
      )
      loss <- res$value
      μ <- res$par
    }
  }
  return(c(loss, τ, μ))
}

updateinfeatures <- function(infeatures, ifit) {
  x <- infeatures
  for (i in ifit) {
    x[i] <- TRUE
  }
  return (x)
}

preparedataSMART <- function(data, param) {

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

  data_standardized <- (data - meanx) / stdx

  return(list(data_standardized = data_standardized, meanx = meanx, stdx = stdx))
}

preparegridsSMART <- function(data, param) {
  τgrid <- gridvectorτ(param$meanlnτ, param$τgridpoints, sharptree = param$sharptree)
  gridmatrixμ <- gridmatrixμ(data, param$μgridpoints)
  μgrid <- gridmatrixμ[[1]]
  dichotomous <- gridmatrixμ[[2]]
  n <- nrow(data)
  p <- ncol(data)

  return(list(τgrid = τgrid, μgrid = μgrid, dichotomous = dichotomous, n = n, p = p))
}



fit_one_tree <- function(r, h, x, infeatures, μgrid, dichotomous, τgrid, param) {
  var_wr <- var(r)
  varϵ <- var_wr * (1 - param$R2p)

  n <- nrow(x)
  p <- ncol(x)
  G0 <- matrix(1, n, 1) # initialize G, the matrix of features
  loss0 <- Inf

  yfit0 <- matrix(0, n, 1)
  ifit <- integer()
  μfit <- numeric()
  τfit <- numeric()
  infeaturesfit <- infeatures
  fi2 <- matrix(0, param$depth, 1)
  βfit <- numeric()

  subsamplesize <- round(n * param$subsamplesharevs)

  if (param$subsamplesharevs == 1) {
    ssi <- 1:n
  } else {
    ssi <- sample(1:n, subsamplesize, replace = FALSE) # sub-sample indexes. Sub-sample no reimmission
  }

  for (depth in 1:param$depth) { #  NB must extend G for this to be viable
    # variable selection
    if (param$subsamplesharevs == 1) {
      outputarray <- loopfeatures(r, h, G0, x, ifit, infeaturesfit, μgrid, dichotomous, τgrid, param, varϵ) # loops over all variables
    } else { # Variable selection using a random sub-set of the sample. All the sample is then used in refinement.
      if (length(h) == 1) {
        outputarray <- loopfeatures(r[ssi], h[ssi], G0[ssi,], x[ssi,], ifit, infeaturesfit, μgrid, dichotomous, τgrid, param, varϵ) # loops over all variables
      } else {
        outputarray <- loopfeatures(r[ssi], h, G0[ssi,], x[ssi,], ifit, infeaturesfit, μgrid, dichotomous, τgrid, param, varϵ) # loops over all variables
      }
    }

    i <- which.min(outputarray[, 1]) # outputarray[,1] is loss (minus log marginal likelihood) vector
    τ0 <- outputarray[i, 2]
    μ0 <- outputarray[i, 3]

    infeaturesfit <- updateinfeatures(inferencesfit, i)

    # refine optimization, after variable selection
    if (param$subsamplesharevs < 1 && param$subsamplefinalbeta == TRUE) {
      if (length(h) == 1) {
        loss <- refineOptim(r[ssi], h[ssi], G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, μ0, dichotomous[i], τ0, param, varϵ)
      } else {
        loss <- refineOptim(r[ssi], h, G0[ssi,], x[ssi, i], infeaturesfit, dichotomous, μ0, dichotomous[i], τ0, param, varϵ)
      }
    } else {
      loss <- refineOptim
    }
  }
}

