# Package
library(deSolve)
library(ggplot2)
library(reshape2)
library(KernSmooth)
library(leaps)
library(snowfall)
library(splines)
library(fda)
library(TruncatedNormal)
library(lpridge)
library(expm)
library(nprobust)
library(Matrix)
library(CollocInfer)
library(magi)

# ==============================================================================
# Data generation
dat_gen <- function(p, J, para, rat, int_upp, omega){
  
  dat_mean <- dat_mean_gen(p, 1000, para, int_upp, omega)[[2]]
  sigma <- sqrt(colSums(dat_mean ^ 2) / 1000) * rat 
  
  dat_mean <- dat_mean_gen(p, J, para, int_upp, omega)
  x <- dat_mean[[1]]
  dat_mean <- matrix(dat_mean[[2]], length(x), p)
  
  dat <- sapply(1:p, function(i){
    dat_mean[,i] + rnorm(J, 0, sigma[i])
  })
  
  colnames(dat) <- colnames(dat_mean)
  return(list(dat = dat, dat_mean = dat_mean, sigma = sigma, para = para, rat = rat, p = p, J = J, x = x))
}

# ==============================================================================
# Estimation
## Bandwidth selection for local polynomial fitting
bandwidth_sel <- function(dat_i, h, sam_grid, ord, k, int_upp){
  mark <- 1:length(sam_grid)
  mark <- mark[(sam_grid >= int_upp / 40) & (sam_grid <= int_upp - int_upp / 40)]
  loss <- sum(sapply(mark, function(j){
    (lpepa(sam_grid[-j], dat_i[-j], bandwidth = h, order = k + 1, deriv = 0, x.out = sam_grid[j])$est - dat_i[j]) ^ 2
  }))
  if(is.finite(loss) == F){
    loss <- 10 ^ 10
  }
  return(loss)
}

## Gradient matching
gra_mat <- function(X, Y, Z, para, weigth){
  opt <- NULL
  init <- 1
  while ((is.null(opt) == T) & (init < 10)) {
    init <- init + 1
    if(init > 2){
      para <- runif(length(para), 0.01, 10)
    }
    
    opt <- sapply(1:20, function(k){
      opt <- optim(par = para, 
                   fn = function(Y, X, Z, para, weigth){
                     sum(sapply(1:length(Y), function(i){
                       sum(((Y[[i]] - X[[i]] %*% con_coef(para, i) + Z[[i]]) * weigth) ^ 2)
                     }))
                   }, 
                   gr = function(Y, X, Z, para, weigth){
                     rowSums(sapply(1:length(Y), function(i){
                       as.numeric(t((Y[[i]] - X[[i]] %*% con_coef(para, i) + Z[[i]]) * weigth) %*%
                                    (X[[i]] * weigth) %*% (jocobi(para, i)) * (-2))
                     }))
                   },
                   Y = Y, X = X, Z = Z,
                   weigth = weigth,
                   method = "L-BFGS-B",
                   lower = rep(0.01, para_length),
                   upper = rep(10, para_length),
                   control = list(maxit = 1000))
      return(c(opt$value, opt$par))
    })
    opt <- opt[, which.min(opt[1,])]
    if(sum(abs(opt[-1]) >= 10) + (sum(abs(opt[-1]) <= 0.01) >= 1)){
      opt <- NULL
    }
  }
  if(is.null(opt) == T){
    opt <- c(Inf, rep(10, para_length))
  }
  
  return(opt)
}

gra_mat_imp <- function(X, Y, Z, para_length, weigth){
  Result <- sapply(1:10, function(k){
    gra_mat(X, Y, Z, runif(para_length, 0.01, 10), weigth)
  })
  Result <- Result[-1,which.min(Result[1,])]
  return(Result)
}

## Integral-equation-based matching
int_mat <- function(p, para, null, XX, Y, Z, cal_term){
  
  para <- optim(par = para, 
                fn = function(para, dim, null, XX, Y, Z, cal_term){
                  val <- sum(sapply(1:dim, function(i){
                    fit <- lm((Y[[i]] - XX[[i]] %*% con_coef(para, i) + Z[[i]]) ~ null[[i]] + 0)
                    return(sum(as.numeric(fit$residuals) ^ 2))
                  }))
                  if(is.finite(val) == F){
                    val <- 10 ^ 10
                  }
                  return(val)
                },
                gr = function(para, dim, null, XX, Y, Z, cal_term){
                  val <- rowSums(sapply(1:dim, function(i){
                    fit <- lm((Y[[i]] - XX[[i]] %*% con_coef(para, i) + Z[[i]]) ~ null[[i]] + 0)
                    return(as.vector(t(as.numeric(fit$residuals)) %*% cal_term[[i]] %*% (jocobi(para, i)) * (-2)))
                  }))
                  return(val)
                },
                dim = p, null = null, XX = XX,
                Y = Y, Z = Z,
                cal_term = cal_term,
                method = "L-BFGS-B",
                lower = rep(0.01, length(para)),
                upper = rep(10, length(para)),
                control = list(maxit = 1000))
  
  return(c(para$value, para$par))
} 

int_con <- function(sam_grid, time_grid, ord){
  if(ord == 1){
    t(sapply(sam_grid, function(j){
      vec <- rep(0, length(time_grid))
      if(j > 0){
        mark <- which(time_grid <= j)
        vec[1:length(mark)] <- 1
      }
      return(vec)
    }))
  }else{
    t(sapply(sam_grid, function(j){
      vec <- rep(0, length(time_grid))
      if(j > 0){
        mark <- which(time_grid <= j)
        vec[1:length(mark)] <- (time_grid[mark] - j) ^ (ord - 1) /  cumprod(1:(ord - 1))[ord - 1]
      }
      return(vec)
    }))
  }
}

int_mat_imp <- function(p, para, null, XX, Y, Z, cal_term){
  Result <- int_mat(p, para, null, XX, Y, Z, cal_term)
  for(ii in 1:10){
    Result <- cbind(Result, int_mat(p, runif(para_length, 0.01, 10), null, XX, Y, Z, cal_term))
  }
  
  Result <- Result[-1,which.min(Result[1,])]
  
  return(Result)
}

prof <- function(int_upp, J, sam_grid, dat, para){
  
  prof_time <- Sys.time()
  bbasis <- create.bspline.basis(
    range = c(0, int_upp),
    nbasis = J,
    norder = 6
  )

  DEfd <- smooth.basis(sam_grid, cbind(dat, rep(0, nrow(dat)), rep(0, nrow(dat))), fdPar(bbasis, 1, 0.1))
  coefs <- DEfd$fd$coefs

  names(par) <- names(para)
  Result <- Profile.LS(
    fn = dynamic_pro,
    data = cbind(dat, rep(NaN, nrow(dat)), rep(NaN, nrow(dat))),
    times = sam_grid,
    pars = para,
    coefs = coefs,
    basisvals = bbasis,
    lambda = 10 ^ 10,
    in.meth = 'nlminb',
    out.meth = 'nls',
    control.out = list(trace = F, printEval = F, warnOnly = T)
  )$par
  
  return(list(Result = Result, prof_time = difftime(Sys.time(), prof_time, units = "secs") ))
}

MAGI <- function(Result, h){
  
  magi_time <- Sys.time()
  sam_grid <- Result[[h]]$sam_grid
  dat <- Result[[h]]$dat
  smo_dat_sam <- Result[[h]]$smo_dat_sam
  para <- Result[[h]]$para
  
  # MAGI
  yinput <- data.frame(sam_grid, dat, rep(NaN, nrow(dat)), rep(NaN, nrow(dat)))
  colnames(yinput)[1] <- c("time")
  
  FNres <- MagiSolver(y = yinput, odeModel = fnmodel,
                      control = list(nstepsHmc = 100, niterHmc = 3000, thetaInit = para,
                                     xInit = cbind(smo_dat_sam[,,1], smo_dat_sam[,,2])))
  return(list(Result = as.vector(summary(FNres)[1,]), magi_time = difftime(Sys.time(), magi_time, units = "secs")))
}

## Parameter estimation for simulation
para_est <- function(dat, ord, para_length, int_upp, omega, para, sam_grid, diff, linear) {
  dat <- as.matrix(dat)
  p <- ncol(dat)
  J <- nrow(dat)
  
  # Pre-smooth
  Init_time <- Sys.time()
  time_grid <- seq(0, int_upp, length.out = 1000)
  smo_dat <- sapply(0:ord, function(k){
    sapply(1:p, function(i){
      h <- dpill(sam_grid, dat[,i])
      h <- optim(h, fn = bandwidth_sel, dat_i = dat[,i], int_upp = int_upp, sam_grid = sam_grid, ord = ord, k = k, method = "Brent",
                 lower = 0, upper = int_upp / 2)$par
      y <- lpepa(sam_grid, dat[,i], bandwidth = h, order = k + 1, deriv = k, x.out = time_grid)$est
      return(y)
    }, simplify = "array")
  }, simplify = "array")
  
  # The fitted values for the sample points
  mark <- sapply(1:length(sam_grid), function(k) which.min(abs(time_grid - sam_grid[k])))
  smo_dat_sam <- smo_dat[mark,,]
  weigth <- (time_grid >= 0) * (time_grid <= int_upp)
  
  Init_time <- difftime(Sys.time(), Init_time, units = "secs")
  
  if(linear == F){
    # Gradient matching
    X <- lapply(1:p, function(i){
      con_state(matrix(smo_dat[,,1:ord], nrow = dim(smo_dat)[1]), i, time_grid) 
    })
    Y <- lapply(1:p, function(i){
      smo_dat[,i,ord + 1]
    })
    if(ord == 1){
      Z <- lapply(1:p, function(i){
        smo_dat[,i,1] * omega[i]
      })
    }else{
      Z <- lapply(1:p, function(i){
        smo_dat[,i,1:ord] %*% omega[i,]
      })
    }
    
    Result_1 <- gra_mat_imp(X, Y, Z, para_length, weigth)
    
    if(sum(abs(Result_1) >= 10) + sum(abs(Result_1) <= 0.01) > 0){
      Result_2 <- rep(10, para_length)
      Result_3 <- rep(10, para_length)
      Result_4 <- rep(10, para_length)
    }else{
      
      # Integral matching
      null <- lapply(1:p, function(i){
        rep(1, length(time_grid)) * weigth
      })
      
      int <- int_con(time_grid, time_grid, 1) * (time_grid[2] - time_grid[1])
      XX <- lapply(1:p, function(i){
        tcrossprod(int, t(X[[i]])) * weigth
      })
      
      Y <- lapply(1:p, function(i){
        smo_dat[,i,ord] * weigth
      })
      
      Z <- lapply(1:p, function(i){
        (cumsum(Z[[i]]) * (time_grid[2] - time_grid[1])) * weigth
      })
      
      cal_term <- lapply(1:p, function(i){
        mat <- null[[i]] %*% solve(t(null[[i]]) %*% null[[i]]) %*% t(null[[i]])
        return(- mat %*% XX[[i]] + XX[[i]])
      })
      
      Result_2 <- int_mat_imp(p, Result_1, null, XX, Y, Z, cal_term)
      
      # Green matching
      if((ord == 1) & (sum(omega != 0) == 0)){
        Result_3 <- Result_2
      }else{
        if(ord == 1){
          null <- lapply(1:p, function(i){
            sapply(time_grid, function(t){
              exp(- t * omega[i])
            }) * weigth
          })
          
          gre <- lapply(1:p, function(i){
            green <- t(sapply(time_grid, function(j){
              vec <- rep(0, length(time_grid))
              if(j > 0){
                mark <- which(time_grid <= j)
                vec[1:length(mark)] <- exp(-(j - time_grid[mark]) * omega[i])
              }
              return(vec)
            }))
            
            return(green)
          })
          
        }else{
          A <- lapply(1:p, function(i){
            rbind(cbind(rep(0, ord - 1), diag(1, ord - 1)), - omega[i,])
          })
          null <- lapply(1:p, function(i){
            t(sapply(time_grid, function(t){
              expm(t * A[[i]], method = "Ward77")[1,]
            }))  * weigth
          })
          
          gre <- lapply(1:p, function(i){
            if(sum(abs(omega[i,])) != 0){
              eigen_A <- eigen(A[[i]])
              eigen_A$svectors <- solve(eigen_A$vectors)
              eigen_A$values <- exp(eigen_A$values)
              
              green <- t(sapply(time_grid, function(j){
                vec <- rep(0, length(time_grid))
                if(j > 0){
                  mark <- (which(time_grid <= j) - 1)
                  vec[length(mark):1] <- sapply(mark, function(t){
                    sum(eigen_A$vectors[1,] * (eigen_A$values ^ (t * (time_grid[2] - time_grid[1]))) * eigen_A$svectors[,2])
                  })
                }
                return(Re(vec))
              }))
            }else{
              green <- t(sapply(time_grid, function(j){
                vec <- rep(0, length(time_grid))
                if(j > 0){
                  mark <- which(time_grid <= j)
                  vec[1:length(mark)] <- (j - time_grid[mark]) ^ (ord - 1)
                }
                return(vec)
              }))
            }
            return(green)
          })
        }
        
        if(diff == T){
          XX <- lapply(1:p, function(i){
            tcrossprod(gre[[i]], t(X[[i]])) * (time_grid[2] - time_grid[1]) * weigth
          })
        }else{
          XX <- lapply(1:p, function(i){
            cbind(tcrossprod(gre[[i]], t(X[[i]][,1:(ncol(X[[i]]) - ord + 1)])) * (time_grid[2] - time_grid[1]) * weigth,
                  sapply(1:(ord-1), function(k){
                    tcrossprod(int_con(time_grid, time_grid, ord - k), t(smo_dat[,i,1])) * (time_grid[2] - time_grid[1]) * weigth
                  }))
          })
        }
        
        Y <- lapply(1:p, function(i){
          smo_dat[,i,1] * weigth
        })
        
        Z <- lapply(1:p, function(i){
          rep(0, length(time_grid))
        })
        
        cal_term <- lapply(1:p, function(i){
          mat <- null[[i]] %*% solve(t(null[[i]]) %*% null[[i]]) %*% t(null[[i]])
          return(- mat %*% XX[[i]] + XX[[i]])
        })
        
        Result_3 <- int_mat_imp(p, Result_1, null, XX, Y, Z, cal_term)
      }
    }
  }else{
    
    gra_time <- Sys.time()
    # Gradient matching
    X <- lapply(1:p, function(i){
      con_state(matrix(smo_dat[,,1:ord], nrow = dim(smo_dat)[1]), i, time_grid) 
    })
    Y <- lapply(1:p, function(i){
      smo_dat[,i,ord + 1]
    })
    if(ord == 1){
      Z <- lapply(1:p, function(i){
        smo_dat[,i,1] * omega[i]
      })
    }else{
      Z <- lapply(1:p, function(i){
        smo_dat[,i,1:ord] %*% omega[i,]
      })
    }
    
    Result_1 <- lapply(1:p, function(i){
      lm(Y[[i]] + Z[[i]] ~ (X[[i]]) + 0)$coefficients
    })
    Result_1 <- as.vector(unlist(Result_1))
    gra_time <- difftime(Sys.time(), gra_time, units = "secs")
    
    # Integral matching
    Inte_time <- Sys.time()
    null <- lapply(1:p, function(i){
      rep(1, length(time_grid)) * weigth
    })
    
    int <- int_con(time_grid, time_grid, 1) * (time_grid[2] - time_grid[1])
    XX <- lapply(1:p, function(i){
      tcrossprod(int, t(X[[i]])) * weigth
    })
    
    Y <- lapply(1:p, function(i){
      smo_dat[,i,ord] * weigth
    })
    
    Z <- lapply(1:p, function(i){
      (cumsum(Z[[i]]) * (time_grid[2] - time_grid[1])) * weigth
    })
    
    Result_2 <- lapply(1:p, function(i){
      lm(Y[[i]] + Z[[i]] ~ (XX[[i]]) + null[[i]] + 0)$coefficients[1:(ncol(XX[[i]]))]
    })
    Result_2 <- as.vector(unlist(Result_2))
    Inte_time <- difftime(Sys.time(), Inte_time, units = "secs")
    
    # Green matching
    Gree_time <- Sys.time()
    if((ord == 1) & (sum(omega != 0) == 0)){
      Result_3 <- Result_2
    }else{
      if(ord == 1){
        null <- lapply(1:p, function(i){
          sapply(time_grid, function(t){
            exp(- t * omega[i])
          }) * weigth
        })
        
        gre <- lapply(1:p, function(i){
          green <- t(sapply(time_grid, function(j){
            vec <- rep(0, length(time_grid))
            if(j > 0){
              mark <- which(time_grid <= j)
              vec[1:length(mark)] <- exp(-(j - time_grid[mark]) * omega[i])
            }
            return(vec)
          }))
          
          return(green)
        })
        
      }else{
        A <- lapply(1:p, function(i){
          rbind(cbind(rep(0, ord - 1), diag(1, ord - 1)), - omega[i,])
        })
        null <- lapply(1:p, function(i){
          t(sapply(time_grid, function(t){
            expm(t * A[[i]], method = "Ward77")[1,]
          }))  * weigth
        })
        
        gre <- lapply(1:p, function(i){
          if(sum(abs(omega[i,])) != 0){
            eigen_A <- eigen(A[[i]])
            eigen_A$svectors <- solve(eigen_A$vectors)
            eigen_A$values <- exp(eigen_A$values)
            
            green <- t(sapply(time_grid, function(j){
              vec <- rep(0, length(time_grid))
              if(j > 0){
                mark <- (which(time_grid <= j) - 1)
                vec[length(mark):1] <- sapply(mark, function(t){
                  sum(eigen_A$vectors[1,] * (eigen_A$values ^ (t * (time_grid[2] - time_grid[1]))) * eigen_A$svectors[,2])
                })
              }
              return(Re(vec))
            }))
          }else{
            green <- t(sapply(time_grid, function(j){
              vec <- rep(0, length(time_grid))
              if(j > 0){
                mark <- which(time_grid <= j)
                vec[1:length(mark)] <- (j - time_grid[mark]) ^ (ord - 1)
              }
              return(vec)
            }))
          }
          return(green)
        })
      }
      
      if(diff == T){
        XX <- lapply(1:p, function(i){
          tcrossprod(gre[[i]], t(X[[i]])) * (time_grid[2] - time_grid[1]) * weigth
        })
      }else{
        XX <- lapply(1:p, function(i){
          cbind(tcrossprod(gre[[i]], t(X[[i]][,1:(ncol(X[[i]]) - ord + 1)])) * (time_grid[2] - time_grid[1]) * weigth,
                - sapply(1:(ord-1), function(k){
                  tcrossprod(int_con(time_grid, time_grid, ord - k), t(smo_dat[,i,1])) * (time_grid[2] - time_grid[1]) * weigth
                }))
        })
      }
      
      Y <- lapply(1:p, function(i){
        smo_dat[,i,1] * weigth
      })
      
      Z <- lapply(1:p, function(i){
        rep(0, length(time_grid))
      })
      
      Result_3 <- lapply(1:p, function(i){
        lm(Y[[i]] + Z[[i]] ~ (XX[[i]]) + null[[i]] + 0)$coefficients[1:ncol(XX[[i]])]
      })
      Result_3 <- as.vector(unlist(Result_3))
    }
    Gree_time <- difftime(Sys.time(), Gree_time, units = "secs")
  }
  
  return(list(Result_1 = Result_1, Result_2 = Result_2, Result_3 = Result_3, 
              smo_dat = smo_dat, time_grid = time_grid, smo_dat_sam = smo_dat_sam,
              gra_time = gra_time + Init_time,
              Inte_time = Inte_time + Init_time,
              Gree_time = Gree_time + Init_time
              ))
}

# ==============================================================================
# Simulation
sim_function <- function(seed, p, J, para, rat, ord, para_length, int_upp, omega, diff, linear){
  set.seed(seed * 100)
  
  dat_sam <- dat_gen(p, J, para, rat, int_upp, omega)
  dat <- dat_sam$dat
  sam_grid <- dat_sam$x
  para_es <- para_est(dat, ord, para_length, int_upp, omega, para, sam_grid, diff, linear)
  
  Result <- list(para_est = para_es, 
                 dat_sam = dat_sam,
                 seed = seed)
  
  return(Result)
}

mse <- function(ture, sim){
  res <- ture - sim
  mse <- mean(sqrt(rowMeans(res ^ 2) / ture ^ 2))
  return(c(mse))
}

deresult <- function(Result){
  
  ture_para <- Result[[1]]$dat_sam$para
  
  fit_para_gra <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_1
  })
  
  fit_para_int <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_2
  })
  
  fit_para_gre <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_3
  })
  
  fit_para_mag <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_4$Result
  })
  
  fit_para_pro <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_5$Result
  })
  
  para_gra <- mse(ture_para, fit_para_gra)
  para_int <- mse(ture_para, fit_para_int)
  para_gre <- mse(ture_para, fit_para_gre)
  para_pro <- mse(ture_para, fit_para_pro)
  para_mag <- mse(ture_para, fit_para_mag)
  
  result <- cbind(para_gra, para_int, para_gre, para_pro, para_mag) * 100
  rownames(result) <- c("rmse")
  colnames(result) <- c("Gradient matching", "Integral matching", "Green's matching", 
                        "Generlized smoothing approach", "MAGI")
  return(result)
}

deresult_time <- function(Result){
  
  fit_para_gra <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$gra_time
  })
  
  fit_para_int <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Inte_time
  })
  
  fit_para_gre <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Gree_time
  })
  
  fit_para_mag <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_4$magi_time
  })
  
  fit_para_pro <- sapply(1:length(Result), function(k) {
    Result[[k]]$para_est$Result_5$prof_time
  })
  
  para_gra <- mean(fit_para_gra)
  para_int <- mean(fit_para_int)
  para_gre <- mean(fit_para_gre)
  para_pro <- mean(fit_para_pro)
  para_mag <- mean(fit_para_mag)
  
  result <- cbind(para_gra, para_int, para_gre, para_pro, para_mag) 
  rownames(result) <- c("Running time")
  colnames(result) <- c("Gradient matching", "Integral matching", "Green's matching", 
                        "Generlized smoothing approach", "MAGI")
  return(result)
}
