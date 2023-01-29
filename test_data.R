pacman::p_load(tidymodels)

data <- diamonds %>%
  select(x,y,z) %>%
  mutate(date = seq(as.Date("1990-01-01"), by = "day", length.out = length(.$x)))


SMARTrecipe <- recipes::recipe(data) %>%
  update_role(x, new_role = "outcome") %>%
  update_role(c(y,z), new_role = "predictor") %>%
  update_role(date, new_role = "date")

x <- SMARTrecipe


library(foreach)
library(doParallel)

remove_parallel_backend <- function(){

  env <- foreach:::.foreachGlobals
  rm(list= ls(name=env), pos=env)

}

devtools::load_all()

cl <- makeCluster(detectCores()-20) # Create a cluster with the number of cores available
registerDoParallel(cl)

tic()
test <- SMARTboost(x ~ z+y,data = diamonds)
toc()
predict(test, diamonds)

stopCluster(cl)

remove_parallel_backend()




m <- matrix(rnorm(9), 3, 3)
foreach(i=1:ncol(m), .combine=c) %do%
  mean(m[,i])

# normalize the rows of a matrix in parallel, with parenthesis used to
# force proper operator precedence
# Need to register a parallel backend before this example will run
# in parallel
suppressWarnings(
foreach(i=1:nrow(m), .combine=rbind) %dopar%
  (m[i,] / mean(m[i,]))
)


d <- data.frame(x=1:10, y=rnorm(10))
s <- foreach(d=iterators::iter(d, by='row'), .combine=rbind) %dopar% d
identical(s, d)












