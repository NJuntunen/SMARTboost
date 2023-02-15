#' @export
SMARTparam <- function(
           loss = "L2",
           coeff = NULL,
           verbose = "Off",
           randomizecv = FALSE,
           ncores = 1,
           sharevalidation = 0.30,
           stderulestop = 0.01,
           stopwhenlossup = FALSE,
           lambda = 0.20,
           depth = 4,
           sigmoid = "sigmoidsqrt",
           meanlntau = 0.0,
           varlntau = 1.0,
           doflntau = 5.0,
           varmu = 4.0,
           dofmu = 10.0,
           subsamplesharevs = 1.0,
           subsamplefinalbeta = FALSE,
           subsampleshare_columns = 1.0,
           mugridpoints = 10,
           taugridpoints = 3,
           refineOptimGrid = FALSE,
           xtolOptim = 0.02,
           optimizevs = FALSE,
           sharptree = FALSE,
           ntrees = 2000,
           R2p = 0.898,
           p0 = 1,
           loglikdivide = 1.0,
           overlap = 0){
  list(loss=loss,
       coeff=coeff,
       verbose=verbose,
       randomizecv=randomizecv,
       ncores=ncores,
       sharevalidation=sharevalidation,
       stderulestop=stderulestop,
       stopwhenlossup=stopwhenlossup,
       lambda=lambda,
       depth=depth,
       sigmoid=sigmoid,
       meanlntau=meanlntau,
       varlntau=varlntau,
       doflntau=doflntau,
       varmu=varmu,
       dofmu=dofmu,
       subsamplesharevs=subsamplesharevs,
       subsamplefinalbeta=subsamplefinalbeta,
       subsampleshare_columns=subsampleshare_columns,
       mugridpoints=mugridpoints,
       taugridpoints=taugridpoints,
       refineOptimGrid=refineOptimGrid,
       xtolOptim=xtolOptim,
       optimizevs=optimizevs,
       sharptree =sharptree,
       ntrees = ntrees,
       R2p = R2p,
       p0 = p0,
       loglikdivide = loglikdivide,
       overlap = overlap)

}


structure(
  list(
    assets = c(1,2)
  )
)










