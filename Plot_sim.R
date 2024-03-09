# Setting
setwd("~") # Set your address

# Function
source("Function.R")

# Load results
rat_list <-  c(0.03, 0.05, 0.07)
sam_list <- c(50, 150, 250)

## Gene regulatory network
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
    
    Result <- deresult(Result)
    Result[3,] <- Result[1,] + Result[2,]
    
    Result[,2] <- Result[,1]
    Result[,1] <- NA
    
    Re_table[[i]][[k]] <- data.frame(Value = c(Result), 
                                     Method = c(rep("  Order 2 gradient matching", 3), rep(" Order 1 gradient matching", 3), rep("Green's matching", 3)),
                                     Measure = rep(c(" RBIAS", "RSD", "  RRMSE"), 3),
                                     Rat = rep(paste0(rat * 100, "%"), 9), 
                                     n = rep(J, 9))
  }
}

dat_plot_1 <- rbind(Re_table[[1]][[1]], Re_table[[1]][[2]], Re_table[[1]][[3]],
                  Re_table[[2]][[1]], Re_table[[2]][[2]], Re_table[[2]][[3]],
                  Re_table[[3]][[1]], Re_table[[3]][[2]], Re_table[[3]][[3]])
dat_plot_1$n <- as.factor(dat_plot_1$n)
dat_plot_1 <- cbind(dat_plot_1, rep("  Gene regulatory network", nrow(dat_plot_1)))
colnames(dat_plot_1)[6] <- "Model"

## Spring-mass system
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
    
    Result <- deresult(Result)
    Result[3,] <- Result[1,] + Result[2,]
    
    Re_table[[i]][[k]] <- data.frame(Value = c(Result), 
                                     Method = c(rep("  Order 2 gradient matching", 3), rep(" Order 1 gradient matching", 3), rep("Green's matching", 3)),
                                     Measure = rep(c(" RBIAS", "RSD", "  RRMSE"), 3),
                                     Rat = rep(paste0(rat * 100, "%"), 9), 
                                     n = rep(J, 9))
  }
}

dat_plot_2 <- rbind(Re_table[[1]][[1]], Re_table[[1]][[2]], Re_table[[1]][[3]],
                  Re_table[[2]][[1]], Re_table[[2]][[2]], Re_table[[2]][[3]],
                  Re_table[[3]][[1]], Re_table[[3]][[2]], Re_table[[3]][[3]])
dat_plot_2$n <- as.factor(dat_plot_2$n)
dat_plot_2 <- cbind(dat_plot_2, rep(" Spring-mass system", nrow(dat_plot_2)))
colnames(dat_plot_2)[6] <- "Model"

## Oscillatory dynamic causal model
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
    
    Result <- deresult(Result)
    Result[3,] <- Result[1,] + Result[2,]
    
    Re_table[[i]][[k]] <- data.frame(Value = c(Result), 
                                     Method = c(rep("  Order 2 gradient matching", 3), rep(" Order 1 gradient matching", 3), rep("Green's matching", 3)),
                                     Measure = rep(c(" RBIAS", "RSD", "  RRMSE"), 3),
                                     Rat = rep(paste0(rat * 100, "%"), 9), 
                                     n = rep(J, 9))
  }
}

dat_plot_3 <- rbind(Re_table[[1]][[1]], Re_table[[1]][[2]], Re_table[[1]][[3]],
                    Re_table[[2]][[1]], Re_table[[2]][[2]], Re_table[[2]][[3]],
                    Re_table[[3]][[1]], Re_table[[3]][[2]], Re_table[[3]][[3]])
dat_plot_3$n <- as.factor(dat_plot_3$n)
dat_plot_3 <- cbind(dat_plot_3, rep("Oscillatory dynamic directional model", nrow(dat_plot_3)))
colnames(dat_plot_3)[6] <- "Model"

## Plot
dat_all_plot <- rbind(dat_plot_1, dat_plot_2, dat_plot_3)
library(ggthemes)
ggplot(dat_all_plot[dat_all_plot$Measure != "RSD",]) + 
  geom_bar(aes(x = n, y = Value, group = Method), alpha = 0.3,
           dat_all_plot[dat_all_plot$Measure == "  RRMSE",],
           stat = "identity",
           position = "dodge",
           size = 1) +
  geom_bar(aes(x = n, y = Value, fill = Method), 
           dat_all_plot[dat_all_plot$Measure == " RBIAS",],
           stat = "identity",
           position = "dodge",
           size = 1) +
  scale_fill_manual(values = c("blue", "orange", "red")) +
  facet_wrap(Model ~  Rat, nrow = 3, scales = "free_y") +
  # geom_hline(yintercept = 0, size  = 1, linetype = 2) + 
  labs(x = "Sample size", y = "RBIAS or RSD (%)",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_tufte() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        text = element_text(size = 20),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  geom_blank(aes(y = y.value), data = data.frame(Model = dat_all_plot[dat_all_plot$Measure != "RSD",]$Model, 
                                                 x.value = rep(rep(c(50, 150, 250), 2), 27),  
                                                 y.value= c(rep(c(0, 40), 27), rep(c(0, 90), 27), rep(c(0, 150), 27)))) 


ggsave(paste0("Plot/para_est.pdf"), width = 14, height = 9, dpi = 300)

