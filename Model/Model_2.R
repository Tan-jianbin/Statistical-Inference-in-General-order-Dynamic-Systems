## Setting
p <- 10
ord <- 2

## Parameter
k <- rep(3, p)

para <- c(k)
para_length <- length(para)
int_upp <- 20
omega <- t(sapply(1:p, function(i) c(0, 0)))
diff <- T
linear <- F

## Model label
mol_lab <- "Newt"

## Data mean generating
con_state <- function(state, i, time_grid){
  if(i == 1){
    matrix(cbind(state[,2] - state[,1] - 2, state[,1] - 2, rep(1, nrow(state))), nrow = nrow(state))
  }else if(i <= 9){
    matrix(cbind(state[,i+1] - state[,i] - 2, state[,i] - state[,i-1] - 2, rep(1, nrow(state))), nrow = nrow(state))
  }else{
    matrix(cbind(state[,i] - state[,i-1] - 2, rep(1, nrow(state)), 5 * sin(time_grid)), nrow = nrow(state))
  }
}

con_coef <- function(para, i){
  if(i <= 9){
    c(para[i + 1], -para[i], 9.8)
  }else{
    c(-para[i], 9.8, 1)
  }
}

jocobi <- function(para, i){
  if(i <= 9){
    gra <- matrix(0, 3, length(para))
    gra[1,i + 1] <- 1
    gra[2,i] <- -1
  }else{
    gra <- matrix(0, 3, length(para))
    gra[1,i] <- -1
  }
  return(gra)
}

## Dynamic equation for simulation
force_func <- function(state, para, p, t){
  sapply(1:p, function(i) con_state(state, i, t) %*% con_coef(para, i))
}

dynamic_equ <- function(t, state, parms, dim, omega) {
  state <- matrix(state, 1)
  ddy <- c(force_func(matrix(state[,1:dim], 1), parms, dim, t)) - rowSums(matrix(state, dim, 2) * omega)
  dy <- state[,(dim+1):(2*dim)]
  return(list(c(dy, ddy)))
}

dat_mean_gen <- function(p, J, para, int_upp, omega){
  x <- seq(0, int_upp, length.out = J)
  initial <- c(seq(2, 20, 2), rep(0, p))
  dat_mean <- ode(y = initial, times = x, func = dynamic_equ, 
                  parms = para, dim = p, omega = omega)
  dat_mean <- dat_mean[,2:(p+1)]
  return(list(x, dat_mean))
}
