## Setting
p <- 50
ord <- 2

## Parameter
para <- c(sapply(1:p, function(k) c(0.5 + seq(-0.1, 0.1, length.out = p)[k],
                                    0.2 * (-1) ^ (k+1), 0.1,
                                    1 + seq(-0.1, 0.1, length.out = p)[k]
)))[-198]
para_length <- length(para)
int_upp <- 20
omega <- t(sapply(1:p, function(i) c(0, 0)))
diff <- F
linear <- T

## Model label
mol_lab <- "brac"

## Data mean generating
con_state <- function(state, i, time_grid){
  if(i <= 49){
    matrix(cbind(-state[,i], state[,i + 1], (time_grid >= 2) * (time_grid <= 3), -state[,i + 50]), nrow = nrow(state))
  }else{
    matrix(cbind(-state[,i], (time_grid >= 2) * (time_grid <= 3), -state[,i + 50]), nrow = nrow(state))
  }
}

con_coef <- function(para, i){
  if(i <= 49){
    c(para[4*i-3], para[4*i -2], para[4*i-1], para[4*i])
  }else{
    c(para[4*i-3], para[4*i-2], para[4*i-1])
  }
}

jocobi <- function(){}

## Dynamic equation for simulation
force_func <- function(state, para, p, t){
  sapply(1:p, function(i) con_state(state, i, t) %*% con_coef(para, i))
}

dynamic_equ <- function(t, state, parms, dim, omega) {
  state <- matrix(state, 1)
  ddy <- c(force_func(matrix(state, 1), parms, dim, t)) - rowSums(matrix(state, dim, 2) * omega)
  dy <- state[,(dim+1):(2*dim)]
  return(list(c(dy, ddy)))
}

dat_mean_gen <- function(p, J, para, int_upp, omega){
  x <- seq(0, int_upp, length.out = J)
  initial <- rep(0, 2 * p)
  dat_mean <- ode(y = initial, times = x, func = dynamic_equ, 
                  parms = para, dim = p, omega = omega)
  dat_mean <- matrix(dat_mean[,2:(p+1)], nrow = nrow(dat_mean))
  return(list(x, dat_mean))
}

