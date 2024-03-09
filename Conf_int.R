# Setting
setwd("~") # Set your address

# Loading code
source("Model/Model_3.R")
source("Function.R")

# Functions
para_est_1 <- function(dat, ord, para_length, int_upp, omega, para, 
                       sam_grid, diff, linear, h, sub_mark) {
  dat <- as.matrix(dat)
  p <- ncol(dat)
  J <- nrow(dat)
  
  # Pre-smooth
  time_grid <- seq(0, int_upp, length.out = 1000)
  
  smo_dat <- array(0, c(length(time_grid), p, ord + 1))
  for(i in sub_mark){
    smo_dat[,i,] <- sapply(0:ord, function(k){
      y <- lpepa(sam_grid, dat[,i], bandwidth = h[i,k+1], order = k + 1, deriv = k, x.out = time_grid)$est
      return(y)
    }, simplify = "array")
    if(i <= p-1){
      smo_dat[,i+1,1] <- lpepa(sam_grid, dat[,i+1], bandwidth = h[i+1,1], order = 1, deriv = 0, x.out = time_grid)$est
    }
  }
  
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
    X <- lapply(sub_mark, function(i){
      con_state(matrix(smo_dat[,,1:ord], nrow = dim(smo_dat)[1]), i, time_grid) 
    })
    Y <- lapply(sub_mark, function(i){
      smo_dat[,i,ord + 1]
    })
    if(ord == 1){
      Z <- lapply(sub_mark, function(i){
        smo_dat[,i,1] * omega[i]
      })
    }else{
      Z <- lapply(sub_mark, function(i){
        smo_dat[,i,1:ord] %*% omega[i,]
      })
    }
    
    Result_1 <- lapply(1, function(i){
      lm(Y[[i]] + Z[[i]] ~ (X[[i]]) + 0)$coefficients
    })
    Result_1 <- as.vector(unlist(Result_1))
    
    # Integral matching
    null <- lapply(sub_mark, function(i){
      rep(1, length(time_grid)) * weigth
    })
    
    int <- int_con(time_grid, time_grid, 1) * (time_grid[2] - time_grid[1])
    XX <- lapply(1, function(i){
      tcrossprod(int, t(X[[i]])) * weigth
    })
    
    Y <- lapply(sub_mark, function(i){
      smo_dat[,i,ord] * weigth
    })
    
    Z <- lapply(1, function(i){
      (cumsum(Z[[i]]) * (time_grid[2] - time_grid[1])) * weigth
    })
    
    Result_2 <- lapply(1, function(i){
      lm(Y[[i]] + Z[[i]] ~ (XX[[i]]) + null[[i]] + 0)$coefficients[1:(ncol(XX[[i]]))]
    })
    Result_2 <- as.vector(unlist(Result_2))
    
    # Green matching
    if((ord == 1) & (sum(omega != 0) == 0)){
      Result_3 <- Result_2
    }else{
      if(ord == 1){
        null <- lapply(sub_mark, function(i){
          sapply(time_grid, function(t){
            exp(- t * omega[i])
          }) * weigth
        })
        
        gre <- lapply(sub_mark, function(i){
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
        A <- lapply(sub_mark, function(i){
          rbind(cbind(rep(0, ord - 1), diag(1, ord - 1)), - omega[i,])
        })
        null <- lapply(1, function(i){
          t(sapply(time_grid, function(t){
            expm(t * A[[i]], method = "Ward77")[1,]
          }))  * weigth
        })
        
        gre <- lapply(1, function(i){
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
        XX <- lapply(1, function(i){
          tcrossprod(gre[[i]], t(X[[i]])) * (time_grid[2] - time_grid[1]) * weigth
        })
      }else{
        XX <- lapply(sub_mark, function(i){
          cbind(tcrossprod(gre[[1]], t(X[[1]][,1:(ncol(X[[1]]) - ord + 1)])) * (time_grid[2] - time_grid[1]) * weigth,
                - sapply(1:(ord-1), function(k){
                  tcrossprod(int_con(time_grid, time_grid, ord - k), t(smo_dat[,i,1])) * (time_grid[2] - time_grid[1]) * weigth
                }))
        })
      }
      
      Y <- lapply(sub_mark, function(i){
        smo_dat[,i,1] * weigth
      })
      
      Z <- lapply(sub_mark, function(i){
        rep(0, length(time_grid))
      })
      
      Result_3 <- lapply(1, function(i){
        lm(Y[[i]] + Z[[i]] ~ (XX[[i]]) + null[[i]] + 0)$coefficients[1:ncol(XX[[i]])]
      })
      Result_3 <- as.vector(unlist(Result_3))
    }
  }
  
  return(list(Result_1 = Result_1, Result_2 = Result_2, Result_3 = Result_3, 
              smo_dat = smo_dat, time_grid = time_grid, smo_dat_sam = smo_dat_sam))
}

Conf <- function(seed, p, J, para, rat, ord, para_length, 
                 int_upp, omega, diff, linear){
  
  set.seed(2022 + seed * 100)
  dat_sam <- dat_gen(p, J, para, rat, int_upp, omega)
  
  # Bootstrap
  sam_grid <- dat_sam$x
  dat <- dat_sam$dat
  
  sub_mark <- 1
  h <- sapply(0:ord, function(k){
    sapply(1:p, function(i){
      if((i - sub_mark >= 0) & (i - sub_mark <= 1)){
        h <- dpill(sam_grid, dat[,i])
        h <- optim(h, fn = bandwidth_sel, dat_i = dat[,i], sam_grid = sam_grid, ord = ord, k = k, method = "Brent",
                   lower = 0, upper = int_upp / 2)$par
      }else{
        h <- 0
      }
      return(h)
    }, simplify = "array")
  }, simplify = "array")
  
  sigma <- sapply(1:p, function(i){
    if((i - sub_mark >= 0) & (i - sub_mark <= 1)){
      B <- t(sapply(1:J, function(x){
        c(1, sam_grid[x])
      }))
      W <- t(sapply(1:J, function(x){
        W <- 3 / 4 * (1 - (sam_grid[x] - sam_grid) ^ 2 /  h[i,1] ^ 2)
        W[W <= 0] <- 0
        W <- diag(W)
        return(t(c(1, sam_grid[x])) %*% solve(t(B) %*% W %*% B) %*% t(B) %*% W)
      }))
      
      sigma <- sqrt(crossprod(W %*% dat[,i] - dat[,i]) / (J - 2 * sum(diag(W)) + sum(W ^ 2)))
    }else{
      sigma <- 0
    }
    return(sigma)
  }, simplify = "array")
  
  Estimated <- para_est_1(dat, ord, para_length, int_upp, omega, para, sam_grid, diff, linear, h, sub_mark)
  smo_dat_sam <- Estimated$smo_dat_sam[,,1]
  Estimated <- cbind(Estimated$Result_1, Estimated$Result_2, Estimated$Result_3)
  
  Boot <- sapply(1:1000, function(ii){
    dat_boot <- sapply(1:p, function(i){
      set.seed(ii * 20200 + seed + i * 101)
      if((i - sub_mark >= 0) & (i - sub_mark <= 1)){
        smo_dat_sam[,i] + rnorm(J, 0, sigma[i])
      }else{
        smo_dat_sam[,i]
      }
    })
    para_es_boot <- para_est_1(dat_boot, ord, para_length, int_upp, omega, para, sam_grid, diff, linear, h, sub_mark)
    Estimated <- cbind(para_es_boot$Result_1, para_es_boot$Result_2, para_es_boot$Result_3)
    return(Estimated)
  }, simplify = "array")
  
  Result <- list(Low = sapply(1:3, function(k) 
                   sapply(1:4, function(l){
                     quantile(Boot[l,k,], pnorm(qnorm((mean(Boot[l,k,] <= Estimated[l,k]))) * 2 + c(qnorm(0.025))))
                   })
                  ),
                 Upper = sapply(1:3, function(k) 
                   sapply(1:4, function(l){
                     quantile(Boot[l,k,], pnorm(qnorm((mean(Boot[l,k,] <= Estimated[l,k]))) * 2 + c(qnorm(0.975))))
                   })
                 ),
                 Boot = Boot)
  
  return(Result)
}

# Computing confidence intervals
J <- 250
rat <- 0.07 ## Change this for obtaining the results with other contaminating levels.
sfInit(parallel = T, cpus = 40)
sfExport("p", "J", "para", "ord", "para_length", "omega", "int_upp", "rat", "diff", "linear")
sfSource("Function.R")
sfExport("con_coef", "con_state", "dat_mean_gen", "force_func", "dynamic_equ",
         "jocobi", "Conf", "para_est_1")
Result_conf <- sfLapply(1:100, Conf,
                   p = p, J = J,
                   rat = rat, para = para, ord = ord, linear = linear,
                   para_length = para_length, int_upp = int_upp,
                   omega = omega, diff = diff)
sfStop()
save(Result_conf, file = paste0("Result/Result_conf_", rat, ".rda"), version = 2)

# Plot
dat_plot <- lapply(1:4, function(kk){
  dat_plot <- sapply(1:100, function(k) Result_conf[[k]]$Low[kk,])
  row.names(dat_plot) <- c("  Order 2 gradient matching", "  Order 1 gradient matching", "Green's matching")
  dat_plot <- melt(dat_plot)
  
  dat_plot_1 <- sapply(1:100, function(k) Result_conf[[k]]$Upper[kk,])
  row.names(dat_plot_1) <- c("  Order 2 gradient matching", "  Order 1 gradient matching", "Green's matching")
  dat_plot_1 <- melt(dat_plot_1)
  dat_plot <- cbind(dat_plot, dat_plot_1$value, rep(para[kk], nrow(dat_plot)))
  colnames(dat_plot) <- c("Method", "x", "Low", "Upp", "Ture")
  dat_plot$x <- as.factor(dat_plot$x)
  dat_plot <- as.data.frame(dat_plot)
  cover <- rep(F, nrow(dat_plot))
  for(i in 1:nrow(dat_plot)){
    if((dat_plot$Low[i] <= dat_plot$Ture[i]) & (dat_plot$Upp[i] >= dat_plot$Ture[i])){
      cover[i] <- T
    }
  }
  lab <- c(1, 4, 2, 3) 
  dat_plot <- data.frame(dat_plot, Cover = cover, Lab = rep(which(lab == kk), nrow(dat_plot)))
  return(dat_plot)
})

dat_plot <- rbind(dat_plot[[1]], dat_plot[[2]], dat_plot[[3]], dat_plot[[4]])

## Covering probability
res <- sapply(1:4, function(k){
  sapply(1:3, function(l){
    mean(dat_plot[dat_plot$Lab == k,]$Cover[seq(l, 300, 3)])
  })
})
res <- res * 100
res <- matrix(as.integer(res), 3, 4)
rownames(res) <- c("  Order 2 gradient matching", "  Order 1 gradient matching", "Green's matching")
xtable::xtable(res) # Table S2

appender <- function(k){
  lab <- c("$\\textit{a}_", "$\\textit{b}_", "$\\textit{c}_", "$\\textit{d}_")
  TeX(paste0(lab, 1, "$"))
}

library(ggthemes)
library(latex2exp)

ggplot(dat_plot) + 
  geom_errorbar(aes(xmin = Low, y = x, xmax = Upp, color = Method), size = 1.5,
               alpha = 0.8, linetype = 1) +
  geom_vline(aes(xintercept = Ture), color = "black", size = 1, linetype = 1, alpha = 0.8) +
  facet_wrap(Method ~ Lab, ncol = 4, scales = "free_x",
             labeller = as_labeller(appender, default = label_parsed)) +
  scale_color_manual(values = c("blue", "orange", "red")) +
  # scale_alpha_discrete(range = c(0.7, 1)) +
  # geom_hline(yintercept = 0, size  = 1, linetype = 2) + 
  labs(x = "Values", y = "Each simulation",
       title = "",
       colour = "", fill = "", linetype = "") +
  # theme_bw(base_family = "Times") +
  theme_tufte() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size = 17),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_blank(aes(x = x.value), data = data.frame(Lab = rep(c(rep(c(1, 1), 100), rep(c(2, 2), 100), 
                                                               rep(c(3, 3), 100), rep(c(4, 4), 100)), 3), 
                                                 x.value= rep(sapply(1:4, function(k) rep(c(min(dat_plot$Low[dat_plot$Lab ==k]), max(dat_plot$Upp[dat_plot$Lab ==k])), 100)), 3))) 

ggsave(paste0("Plot/CI_", rat, ".pdf"), width = 11.4, height = 11.4, dpi = 300)

