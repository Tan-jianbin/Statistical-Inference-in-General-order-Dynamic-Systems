# Setting
setwd("~") # Set your address

rat_list <-  c(0.03, 0.05, 0.07) # Contamination level
sam_list <- c(50, 150, 250) # Sample size

# ==============================================================================
for(J in sam_list){
  for(rat in rat_list){
    source("Model/Handwriting_sim.R")
    source("Simulation_hand.R")
  }
}

## MAGI
for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "hawr", ".rda"))
    
    Est_magi <- lapply(1:length(Result), function(h){
      list(sam_grid = Result[[h]]$dat_sam$x,
           dat = Result[[h]]$dat_sam$dat,
           smo_dat_sam = Result[[h]]$para_est$smo_dat_sam,
           para = Result[[h]]$para_est$Result_3)
    })
    
    sfInit(parallel = T, cpus = 40)
    sfExport("Est_magi")
    sfSource("Function_hand.R")
    sfSource("Model/Handwriting_sim.R")
    Result_magi <- sfLapply(1:100, MAGI, Result = Est_magi)
    sfStop()
    
    for(h in 1:length(Result)){
      Result[[h]]$para_est$Result_4 <- Result_magi[[h]]
    }
    save(Result, file = paste0("Result/Result_", J, "_", rat, "_", mol_lab, ".rda"), version = 2)
  }
}

## Generalized smoothing approach
for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "hawr", ".rda"))
    for(h in 1:length(Result)){
      Result[[h]]$para_est$Result_5 <- prof(int_upp, J, Result[[h]]$dat_sam$x, Result[[h]]$dat_sam$dat, Result[[h]]$para_est$Result_3)
    }
    save(Result, file = paste0("Result/Result_", J, "_", rat, "_", mol_lab, ".rda"), version = 2)
  }
}

## Results
### Estimated Result
Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "hawr", ".rda"))
    Re_table[[i]][[k]] <- deresult(Result)
  }
}

re_latex <- matrix(0, 15, 3)
for(i in 1:3){
  re_latex[(5*(i-1) + 1):(5*(i-1) + 5),] <- cbind(t(Re_table[[i]][[1]]), t(Re_table[[i]][[2]]), t(Re_table[[i]][[3]]))
}
re_latex <- round(re_latex, 2)
re_latex <- cbind(rep(c("Gradient matching", "Integral matching", "Green’s matching", "Generalized smoothing approach", "MAGI"), 3), re_latex)
rownames(re_latex) <- rep(NULL, 15)
xtable::xtable(re_latex) 

### Running time
Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "hawr", ".rda"))
    Re_table[[i]][[k]] <- deresult_time(Result)
  }
}

re_latex <- matrix(0, 15, 3)
for(i in 1:3){
  re_latex[(5*(i-1) + 1):(5*(i-1) + 5),] <- cbind(t(Re_table[[i]][[1]]), t(Re_table[[i]][[2]]), t(Re_table[[i]][[3]]))
}
re_latex <- round(re_latex, 2)
re_latex <- cbind(rep(c("Gradient matching", "Integral matching", "Green’s matching", "Generalized smoothing approach", "MAGI"), 3), re_latex)
rownames(re_latex) <- rep(NULL, 15)
xtable::xtable(re_latex) 
