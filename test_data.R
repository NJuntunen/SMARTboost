pacman::p_load(tidymodels, profvis, nloptr, tictoc, purrr, SMARTboost)

data <- diamonds %>%
  select(x,y,z) %>%
  mutate(date = seq(as.Date("1990-01-01"), by = "day", length.out = length(.$x)),
         across(where(is.numeric), ~.*100),
         across(where(is.numeric), ~as.integer(.)))


SMARTrecipe <- recipes::recipe(data) %>%
  update_role(x, new_role = "outcome") %>%
  update_role(c(y,z), new_role = "predictor") %>%
  update_role(date, new_role = "date")

x <- SMARTrecipe

tic()
test <- SMARTboost_fit(x ~ z+y,data = diamonds, ntrees = 10, optimizevs = FALSE, ncores = 4)
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




add_boost_tree_SMARTboost()

model_spec <- boost_tree(trees = 10) %>%
  set_mode("regression") %>%
  set_engine("SMARTboost")







