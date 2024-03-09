# Setting
address <- "~"
setwd(address) # Set your address

rat_list <-  c(0.03, 0.05, 0.07) # Contamination level
sam_list <- c(50, 150, 250) # Sample size
  
# ==============================================================================
## Gene regulatory network
for(J in sam_list){
    for(rat in rat_list){
      source("Model/Model_1.R")
      source("Simulation.R")
    }
}

Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "Gene", ".rda"))
    Re_table[[i]][[k]] <- deresult(Result)[, c(1, 3)]
  }
}

re_latex <- matrix(0, 6, 3)
for(i in 1:3){
  re_latex[(2*(i-1) + 1):(2*(i-1) + 2),] <- cbind(t(Re_table[[i]][[1]])[,3], t(Re_table[[i]][[2]])[,3], t(Re_table[[i]][[3]])[,3])
}
re_latex <- round(re_latex, 2)
re_latex <- cbind(rep(c("Gradient matching", "Green’s matching"), 3), re_latex)
rownames(re_latex) <- c(" ", "  ", "   ", "    ", "     ", "      ")
xtable::xtable(re_latex) #  Results in Table 1

# ==============================================================================
## Spring-mass system
for(J in sam_list){
  for(rat in rat_list){
    source("Model/Model_2.R")
    source("Simulation.R")
  }
}

Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "Newt", ".rda"))
    Re_table[[i]][[k]] <- deresult(Result)
  }
}

re_latex <- matrix(0, 9, 3)
for(i in 1:3){
  re_latex[(3*(i-1) + 1):(3*(i-1) + 3),] <- cbind(t(Re_table[[i]][[1]])[,3], t(Re_table[[i]][[2]])[,3], t(Re_table[[i]][[3]])[,3])
}
re_latex <- round(re_latex, 2)
re_latex <- cbind(rep(c("Order 2 Gradient matching", "Order 1 Gradient matching", "Green’s matching"), 3), re_latex)
rownames(re_latex) <- c(" ", "  ", "   ", "    ", "     ", "      ", "       ", "        ", "         ")
xtable::xtable(re_latex) #  Results in Table 1

# ==============================================================================
## Oscillatory dynamic causal model
for(J in sam_list){
  for(rat in rat_list){
    source("Model/Model_3.R")
    source("Simulation.R")
  }
}

Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    J <- sam_list[i]
    rat <- rat_list[k]
    load(paste0("Result/Result_", J, "_", rat, "_", "brac", ".rda"))
    Re_table[[i]][[k]] <- deresult(Result)
  }
}

re_latex <- matrix(0, 9, 3)
for(i in 1:3){
  re_latex[(3*(i-1) + 1):(3*(i-1) + 3),] <- cbind(t(Re_table[[i]][[1]])[,3], t(Re_table[[i]][[2]])[,3], t(Re_table[[i]][[3]])[,3])
}
re_latex <- round(re_latex, 2)
re_latex <- cbind(rep(c("Order 2 gradient matching", "Order 1 gradient matching", "Green’s matching"), 3), re_latex)
rownames(re_latex) <- c(" ", "  ", "   ", "    ", "     ", "      ", "       ", "        ", "         ")
xtable::xtable(re_latex) #  Results in Table 1



