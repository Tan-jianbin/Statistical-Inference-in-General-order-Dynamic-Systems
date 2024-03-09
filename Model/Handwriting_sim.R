# Handwriting
## Setting
p <- 2
ord <- 2

## Parameter
para <- c(-0.0015, c(seq(-0.0005, 0.0005, length.out = 6), -0.0015, seq(-0.002, 0.001, length.out = 6)))
para_length <- length(para)
int_upp <- 600
omega <- t(sapply(1:p, function(i) c(0, 0)))
diff <- T
linear <- T

## Model label
mol_lab <- "hawr"

## Data mean generating
con_state <- function(state, i, time_grid){
  if(i == 1){
    matrix(cbind(state[,1], matrix(sapply(1:6, function(k) (time_grid >= 600 / 6 * (k - 1)) * (time_grid <= 600 / 6 * k)), nrow = nrow(state))), nrow = nrow(state))
  }else{
    matrix(cbind(state[,2], matrix(sapply(1:6, function(k) (time_grid >= 600 / 6 * (k - 1)) * (time_grid <= 600 / 6 * k)), nrow = nrow(state))), nrow = nrow(state))
  }
}

con_coef <- function(para, i){
  if(i == 1){
    para[1:7]
  }else{
    para[8:14]
  }
}

## Dynamic equation for simulation
force_func <- function(state, para, p, t){
  sapply(1:p, function(i) con_state(state, i, t) %*% con_coef(para, i))
}

dynamic_equ <- function(t, state, parms, dim, omega) {
  state <- matrix(state, 1)
  ddy <- c(force_func(matrix(state, 1), parms, dim, t)) - rowSums(matrix(state, dim, 2) * omega)
  dy <- state[,3:4]
  return(list(c(dy, ddy)))
}

dat_mean_gen <- function(p, J, para, int_upp, omega){
  x <- seq(0, int_upp, length.out = J)
  initial <- rep(0, 4)
  dat_mean <- ode(y = initial, times = x, func = dynamic_equ, 
                  parms = para, dim = p, omega = omega)
  dat_mean <- matrix(dat_mean[,2:(p+1)], nrow = nrow(dat_mean))
  return(list(x, dat_mean))
}

# Generalized smoothing approach
fnn <- function(times, x, p, more) {
  r <- x
  r[,1] <- x[,3]
  r[,2] <- x[,4]
  r[,3] <- con_state(x, 1, times) %*% con_coef(p, 1)
  r[,4] <- con_state(x, 2, times) %*% con_coef(p, 2)
  return(r)
}

dfdx <- function(times, x, p, more){
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  resultDx[, 1, 3] <- 1
  resultDx[, 2, 4] <- 1
  resultDx[, 3, 1] <- con_coef(p, 1)[1]
  resultDx[, 4, 2] <- con_coef(p, 2)[1]
  return(resultDx)
}

dfdp <- function(times, x, p, more){
  resultDtheta <- array(0, c(nrow(x), ncol(x), length(p)))
  
  resultDtheta[, 3, 1] <- x[, 1]
  for(k in 1:6){
    resultDtheta[, 3, k + 1] <- (times >= 600 / 6 * (k - 1)) * (times <= 600 / 6 * k)
  }
  
  resultDtheta[, 4, 8] <- x[, 2]
  for(k in 1:6){
    resultDtheta[, 4, k + 8] <- (times >= 600 / 6 * (k - 1)) * (times <= 600 / 6 * k)
  }
  return(resultDtheta)
}

d2fdx2 <- function(times, x, p, more){
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x), ncol(x)))
  return(resultDx)
}

d2fdxdp <- function(times, x, p, more){
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x), length(p)))
  resultDx[, 3, 1, 1] <- 1
  resultDx[, 4, 2, 8] <- 1
  return(resultDx)
}

dynamic_pro <- list(
  fnn = fnn,
  dfdx  = dfdx,
  dfdp = dfdp, 
  d2fdx2 = d2fdx2, 
  d2fdxdp = d2fdxdp
)

# MAGI
fnmodelODE <- function(theta, x, tvec){
  
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] <- x[, 3]
  result[, 2] <- x[, 4]
  result[, 3] <- con_state(x, 1, tvec) %*% con_coef(theta, 1)
  result[, 4] <- con_state(x, 2, tvec) %*% con_coef(theta, 2)
  return(result)
}

fnmodelDx <- function(theta, x, tvec){
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  resultDx[, 3, 1] <- 1
  resultDx[, 4, 2] <- 1
  resultDx[, 1, 3] <- con_coef(theta, 1)[1]
  resultDx[, 2, 4] <- con_coef(theta, 2)[1]
  return(resultDx)
}

fnmodelDtheta <- function(theta, x, tvec){
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  resultDtheta[, 1, 3] <- x[, 1]
  for(k in 1:6){
    resultDtheta[, k + 1, 3] <- (tvec >= 600 / 6 * (k - 1)) * (tvec <= 600 / 6 * k)
  }
  
  resultDtheta[, 8, 4] <- x[, 2]
  for(k in 1:6){
    resultDtheta[, k + 8, 4] <- (tvec >= 600 / 6 * (k - 1)) * (tvec <= 600 / 6 * k)
  }
  return(resultDtheta)
}

fnmodel <- list(
  fOde = fnmodelODE,
  fOdeDx = fnmodelDx,
  fOdeDtheta = fnmodelDtheta,
  thetaLowerBound = rep(-Inf, para_length),
  thetaUpperBound = rep(Inf, para_length)
)