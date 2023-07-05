logit <- function(p) {
  
  x <- log(p/(1-p))
  
}

expit <- function(p) {
  
  1/(1+exp(-p))
  
}


logit_jacobian <- function(x) {

mat <- diag(1/(x-x^2))  
  
}




