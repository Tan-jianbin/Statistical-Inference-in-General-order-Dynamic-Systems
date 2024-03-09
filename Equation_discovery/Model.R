# Loading functions
source("Function.R")

# Setting
p <- 1
ord <- 2

# Parameter
int_upp <- 20
para <- -1
omega <- t(sapply(1:p, function(i) c(0, 0)))

# Data mean generating
con_state <- function(state, time_grid){
  matrix(cbind(sin(state[,1]), cos(state[,1]),
               rep(1, nrow(state)),
               state[,1], state[,1] ^ 2, state[,1] ^ 3, state[,1] ^ 4
               ), nrow = nrow(state))
}

con_coef <- function(para){
  par <- rep(0, 7)
  par[1] <- para
  return(par)
}

# Dynamic equation for simulation
force_func <- function(state, para, p, t){
  sapply(1:p, function(i) con_state(state, t) %*% con_coef(para))
}

dynamic_equ <- function(t, state, parms, dim, omega) {
  state <- matrix(state, 1)
  ddy <- c(force_func(matrix(state, 1), parms, dim, t))
  dy <- state[,(dim+1):(2*dim)]
  return(list(c(dy, ddy)))
}

dat_mean_gen <- function(p, J, para, int_upp, omega){
  x <- seq(0, int_upp, length.out = J)
  initial <- runif(2, -0.5, 0.5)
  dat_mean <- ode(y = initial, times = x, func = dynamic_equ, 
                  parms = para, dim = p, omega = omega)
  dat_mean <- matrix(dat_mean[,2:(p+1)], nrow = nrow(dat_mean))
  return(list(x, dat_mean))
}

## Parameter estimation for simulation
para_est <- function(dat, ord, int_upp, omega, para, sam_grid) {
  dat <- as.matrix(dat)
  p <- ncol(dat)
  J <- nrow(dat)
  
  # Pre-smooth
  time_grid <- seq(0, int_upp, length.out = 1000)
  
  smo_dat <- sapply(0:ord, function(k){
    sapply(1:p, function(i){
      h <- dpill(sam_grid, dat[,i])
      h <- optim(h, fn = bandwidth_sel, dat_i = dat[,i], sam_grid = sam_grid, ord = ord, k = k, method = "Brent",
                 lower = 0, upper = int_upp / 2)$par
      y <- lpepa(sam_grid, dat[,i], bandwidth = h, order = k + 1, deriv = k, x.out = time_grid)$est
      return(y)
    }, simplify = "array")
  }, simplify = "array")
  
  # The fitted values for the sample points
  mark <- sapply(1:length(sam_grid), function(k) which.min(abs(time_grid - sam_grid[k])))
  smo_dat_sam <- smo_dat[mark,,]
  weigth <- (time_grid >= 0.5) * (time_grid <= int_upp - 0.5)
  
  # Gradient matching for equation discovery
  X <- con_state(matrix(smo_dat[,,1:ord], nrow = dim(smo_dat)[1]), time_grid) * weigth
  Y <- smo_dat[,1,ord + 1] * weigth
  
  tun <- cv.glmnet(x = X, y = Y, standardize = F, intercept = F)
  Result_1 <- glmnet(x = X, y = Y, standardize = F, intercept = F, lambda = tun$lambda.min)
  
  # Green's matching for equation discovery
  A <- lapply(1:p, function(i){
    rbind(cbind(rep(0, ord - 1), diag(1, ord - 1)), - omega[i,])
  })
  
  null <- lapply(1:p, function(i){
    t(sapply(time_grid, function(t){
      expm(t * A[[i]], method = "Ward77")[1,]
    }))  * weigth
  })
  
  gre <- lapply(1:p, function(i){
    green <- t(sapply(time_grid, function(j){
      vec <- rep(0, length(time_grid))
      if(j > 0){
        mark <- which(time_grid <= j)
        vec[1:length(mark)] <- (j - time_grid[mark]) ^ (ord - 1)
      }
      return(vec)
    }))
    return(green)
  })
  
  X <- con_state(matrix(smo_dat[,,1:ord], nrow = dim(smo_dat)[1]), time_grid)
  XX <- tcrossprod(gre[[1]], t(X)) * (time_grid[2] - time_grid[1]) * weigth
  XX <- cbind(XX, null[[1]])
  
  Y <- smo_dat[,1,1] * weigth
  
  tun <- cv.glmnet(x = XX, y = Y, standardize = F, intercept = F, penalty.factor = c(rep(1, ncol(X)), 0 ,0))
  Result_2 <- glmnet(x = XX, y = Y, standardize = F, intercept = F, lambda = tun$lambda.min,
                     penalty.factor = c(rep(1, ncol(X)), 0 ,0))
  
  return(list(Result_1 = Result_1, Result_2 = Result_2, 
              smo_dat = smo_dat, time_grid = time_grid, smo_dat_sam = smo_dat_sam))
}

# ==============================================================================
# Simulation
sim_function <- function(seed, p, J, para, rat, ord, int_upp, omega){
  set.seed(seed * 100)
  dat_sam <- dat_gen(p, J, para, rat, int_upp, omega)
  para_es <- para_est(dat = dat_sam$dat, ord, int_upp, omega, para, sam_grid = dat_sam$x)
  
  Result <- list(para = para_es, 
                 dat_sam = dat_sam,
                 seed = seed)
  
  return(Result)
}

