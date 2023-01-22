
initialize_gamma <- function(data, param){

  sum(data, na.rm = T)/length(data)

}

evaluate_pseudoresid <- function(data, gammafit){

  r <- data - gammafit
  return(list(r = r, h = rep(1, length(r))))

}
