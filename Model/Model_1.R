## Setting
p <- 50
ord <- 1 

## Parameter
para <- c(sapply(1:p, function(k){
  c(0.5, 1, 0.1)
}))
 
para_length <- length(para)
int_upp <- 20
omega <- rep(0, p)
diff <- T
linear <- T

## Model label
mol_lab <- "Gene"

## Data mean generating
con_state <- function(state, i, time_grid){
  mark <- c(i, (i + 5) %% 50, (i - 3) %% 50)
  mark[mark == 0] <- 50
  return(matrix(cbind(sin(1  / 2 * state[,mark[1]]), - cos(1 / 2 * state[,mark[2]]), state[,mark[3]]), nrow(state)))
}

con_coef <- function(para, i){
  c(para[3*i-2], para[3*i-1], para[3*i])
}

jocobi <- function(para, i){}

## Dynamic equation for simulation
force_func <- function(state, para, p, t){
  sapply(1:p, function(i) con_state(state, i, t) %*% con_coef(para, i))
}

dynamic_equ <- function(t, state, parms, dim, omega) {
  state <- matrix(state, 1)
  dy <- c(force_func(state, parms, dim, t))
  return(list(dy))
}

dat_mean_gen <- function(p, J, para, int_upp, omega){
  x <- seq(0, int_upp, length.out = J)
  initial <- seq(5, 1, length.out = p)
  dat_mean <- ode(y = initial, times = x, func = dynamic_equ, parms = para, dim = p, omega = omega)
  dat_mean <- dat_mean[,2:(p+1)]
  return(list(x, dat_mean))
}

