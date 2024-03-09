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
library(glmnet)


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
bandwidth_sel <- function(dat_i, h, sam_grid, ord, k){
  mark <- 1:length(sam_grid)
  mark <- mark[(sam_grid >= 0.5) & (sam_grid <= 19.5)]
  loss <- sum(sapply(mark, function(j){
    (lpepa(sam_grid[-j], dat_i[-j], bandwidth = h, order = k + 1, deriv = 0, x.out = sam_grid[j])$est - dat_i[j]) ^ 2
    # (lpridge(sam_grid[-j], dat_i[-j], bandwidth = h, order = k + 1, deriv = 0, x.out = sam_grid[j], ridge = 0, weight = "bi")$est - dat_i[j]) ^ 2
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

## Parameter estimation for simulation
para_est <- function(dat, ord, para_length, int_upp, omega, para, sam_grid, diff, linear) {
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
      # y <- lpridge(sam_grid, dat[,i], bandwidth = h, order = k + 1, deriv = k, x.out = time_grid, ridge = 0, weight = "bi")$est
      y <- lpepa(sam_grid, dat[,i], bandwidth = h, order = k + 1, deriv = k, x.out = time_grid)$est
      return(y)
    }, simplify = "array")
  }, simplify = "array")
  
  # The fitted values for the sample points
  mark <- sapply(1:length(sam_grid), function(k) which.min(abs(time_grid - sam_grid[k])))
  smo_dat_sam <- smo_dat[mark,,]
  weigth <- (time_grid >= 0.5) * (time_grid <= int_upp - 0.5)
  
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
    
    Result_2 <- lapply(1:p, function(i){
      lm(Y[[i]] + Z[[i]] ~ (XX[[i]]) + null[[i]] + 0)$coefficients[1:(ncol(XX[[i]]))]
    })
    Result_2 <- as.vector(unlist(Result_2))
    
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
  }
    
  return(list(Result_1 = Result_1, Result_2 = Result_2, Result_3 = Result_3, 
              smo_dat = smo_dat, time_grid = time_grid, smo_dat_sam = smo_dat_sam))
}

# ==============================================================================
# Simulation
sim_function <- function(seed, p, J, para, rat, ord, para_length, int_upp, omega, diff, linear){
  set.seed(seed * 100)
  para_es <- NULL
  while(is.null(para_es) == T){
    try({
      dat_sam <- dat_gen(p, J, para, rat, int_upp, omega)
      para_es <- para_est(dat_sam$dat, ord, para_length, int_upp, omega, para, dat_sam$x, diff, linear)
      if(sum(abs(para_es$Result_1) == 10) + sum(abs(para_es$Result_1) == 0.01) + 
         sum(abs(para_es$Result_2) == 10) + sum(abs(para_es$Result_2) == 0.01) +
         sum(abs(para_es$Result_3) == 10) + sum(abs(para_es$Result_3) == 0.01) > 0){
        para_es <- NULL
      }
    }, silent = T)
  }
  
  Result <- list(para = para_es, 
                 dat_sam = dat_sam,
                 seed = seed)
  
  return(Result)
}

mse <- function(ture, sim){
  res <- ture - sim
  bias <- mean(abs(rowMeans(res) / ture))
  sd <- mean(sqrt(rowMeans((sim - rowMeans(sim)) ^ 2) / ture ^ 2))
  mse <- mean(sqrt(rowMeans(res ^ 2) / ture ^ 2))
  return(c(bias, sd, mse))
}

deresult <- function(Result){
  
  ture_para <- Result[[1]]$dat_sam$para
  
  fit_para_gra <- sapply(1:length(Result), function(k) {
    Result[[k]]$para$Result_1
  })
  
  fit_para_int <- sapply(1:length(Result), function(k) {
    Result[[k]]$para$Result_2
  })
  
  fit_para_gre <- sapply(1:length(Result), function(k) {
    Result[[k]]$para$Result_3
  })
  
  para_gra <- mse(ture_para, fit_para_gra)
  para_int <- mse(ture_para, fit_para_int)
  para_gre <- mse(ture_para, fit_para_gre)
  
  result <- cbind(para_gra, para_int, para_gre) * 100
  rownames(result) <- c("rbias", "rsd", "rmse")
  colnames(result) <- c("Gradient matching", "Integral matching", "Green's matching")
  return(result)
}
