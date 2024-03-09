# Setting
setwd("~") # Set your address

# Function
source("Function.R")

## Oscillatory dynamic causal model
source("Model/Model_3.R")
J <- 250
rat <- 0.07 ## Change this for obtaining the results with other contaminating levels.
load(paste0("Result/Result_", J, "_", rat, "_", "Brac", ".rda"))

kk <- 1
x <- seq(0, int_upp, length.out = J)
Result <- sapply(1:length(Result), function(ii){
  initial <- c(Result[[ii]]$para$smo_dat[1,,1:ord])
  fit <- array(0, c(3, J, 3))
  
  fit_1 <- ode(y = initial, times = x, func = dynamic_equ, parms = Result[[ii]]$para$Result_1,
               dim = p, omega = omega)[,-1]
  fit[1,,] <- cbind(fit_1[,c(kk, kk+p)],con_state(fit_1, kk, x) %*% con_coef(Result[[ii]]$para$Result_1, kk))
  
  fit_2 <- ode(y = initial, times = x, func = dynamic_equ, parms = Result[[ii]]$para$Result_2,
               dim = p, omega = omega)[,-1]
  fit[2,,] <- cbind(fit_2[,c(kk, kk+p)],con_state(fit_2, kk, x) %*% con_coef(Result[[ii]]$para$Result_2, kk))
  
  fit_3 <- ode(y = initial, times = x, func = dynamic_equ, parms = Result[[ii]]$para$Result_3,
               dim = p, omega = omega)[,-1]
  fit[3,,] <- cbind(fit_3[,c(kk, kk+p)],con_state(fit_3, kk, x) %*% con_coef(Result[[ii]]$para$Result_3, kk))
  
  return(fit)
}, simplify = "array")

fit_1 <- ode(y = rep(0, 2 * p), times = x, func = dynamic_equ, parms = para,
             dim = p, omega = omega)[,-1]
Ture <- cbind(fit_1[,c(kk, kk+p)],con_state(fit_1, kk, x) %*% con_coef(para, kk))

Re_table <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    return(data.frame(fit = c(Result[i,,k,]),
                      x = rep(x, 100),
                      Lab = c(sapply(1:100, function(k) rep(k, 250))),
                      Method = rep(c("  Order 2 gradient matching", " Order 1 gradient matching", "Green's matching")[i], 250 * 100),
                      Diff = rep(c("  Zero derivative", " First derivative", "Second derivative")[k], 250 * 100),
                      Ture = rep(Ture[,k], 100)))
  })
})

## Plot
dat_all_plot <- rbind(Re_table[[1]][[1]], Re_table[[1]][[2]], Re_table[[1]][[3]],
                      Re_table[[2]][[1]], Re_table[[2]][[2]], Re_table[[2]][[3]],
                      Re_table[[3]][[1]], Re_table[[3]][[2]], Re_table[[3]][[3]])

library(ggthemes)
ggplot(dat_all_plot) + 
  geom_line(aes(x = x, y = fit, color = Method, group = Lab),
            alpha = 0.6,
            size = 0.3) +
  geom_line(aes(x = x, y = Ture), 
            size = 0.5) +
  scale_color_manual(values = c("blue", "orange", "red")) +
  facet_wrap(Method ~ Diff, scales = "free_y") +
  # geom_hline(yintercept = 0, size  = 1, linetype = 2) + 
  labs(x = "t", y = "Values",
       title = "",
       colour = "", fill = "", linetype = "") +
  # theme_bw(base_family = "Times") +
  ylim(c(-0.1, 0.1)) +
  theme_tufte() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        text = element_text(size = 20),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) 

ggsave(paste0("Plot/CU_", rat, ".pdf"), width = 11.4, height = 11.4, dpi = 300)

