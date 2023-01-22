SMARTtree <-  function(i, μ, τ, β, fi2) {
  structure(list(i = i,
                 μ = μ,
                 τ = τ,
                 β = β,
                 fi2 = fi2),
            class = c("SMARTtree", "list"))
}

SMARTboostTrees <- function(param, gamma0, n, p, meanx, stdx) {
  structure(list(param = param,
                 gamma0 = gamma0,
                 trees = list(),
                 infeatures = rep(FALSE, p),
                 fi2 = rep(0, p),
                 meanx = meanx,
                 stdx = stdx,
                 gammafit = rep(gamma0, n),
                 R2simul = numeric()),
            class = c("SMARTboostTrees", "list"))
}

updateSMARTtrees <- function(SMARTtrees, Gβ, tree, rh, iter) {

  n <- length(Gβ)
  depth <- length(tree$i)

  SMARTtrees$gammafit <- SMARTtrees$gammafit + SMARTtrees$param$lambda * Gβ
  SMARTtrees$infeatures <- updateinfeatures(SMARTtrees$infeatures, tree$i)
  SMARTtrees$trees <- append(SMARTtrees$trees, list(tree))

  for (d in 1:depth) {
    SMARTtrees$fi2[tree$i[d]] <- SMARTtrees$fi2[tree$i[d]] + tree$fi2[d]
  }

  if (iter == 1 && SMARTtrees$param$R2p == as.numeric(T(0.898))) {
    sqrth <- sqrt(rh$h)
    R2tree <- var(Gβ * sqrth) / var(rh$r / sqrth)
    SMARTtrees$param$R2p <- R2tree
  }

  SMARTtrees
}
