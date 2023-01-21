SMARTparam <- function(
           loss = "L2",
           coeff = NULL,
           verbose = "Off",
           randomizecv = FALSE, 
           nfold = 5,
           sharevalidation = 0.30,
           stderulestop = 0.01,
           stopwhenlossup = FALSE,
           lambda = 0.20,
           depth = 4, 
           sigmoid = "sigmoidsqrt",
           meanlnτ = 0.0,
           varlnτ = 1.0,
           doflnτ = 5.0,
           varμ = 4.0,
           dofμ = 10.0,
           subsamplesharevs = 1.0, 
           subsamplefinalbeta = FALSE,
           subsampleshare_columns = 1.0,
           μgridpoints = 10,
           τgridpoints = 3,
           refineOptimGrid = FALSE, 
           xtolOptim = 0.02,
           optimizevs = FALSE){
  list(loss=loss,
       coeff=coeff,
       verbose=verbose,
       randomizecv=randomizecv,
       nfold=nfold,
       sharevalidation=sharevalidation,
       stderulestop=stderulestop,
       stopwhenlossup=stopwhenlossup,
       lambda=lambda,
       depth=depth,
       sigmoid=sigmoid,
       meanlnτ=meanlnτ,
       varlnτ=varlnτ,
       doflnτ=doflnτ,
       varμ=varμ,
       dofμ=dofμ,
       subsamplesharevs=subsamplesharevs,
       subsamplefinalbeta=subsamplefinalbeta,
       subsampleshare_columns=subsampleshare_columns,
       μgridpoints=μgridpoints,
       τgridpoints=τgridpoints,
       refineOptimGrid=refineOptimGrid,
       xtolOptim=xtolOptim,
       optimizevs=optimizevs,
       sharptree=sharptree,
       ntrees=ntrees,
       R2p=R2p)
  
}














