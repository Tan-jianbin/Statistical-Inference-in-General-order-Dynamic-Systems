# Set address
setwd("~")

# Function
source("Function_data_illu.R")

# Load data
load("data/dat.rda")
source("Model/Handwriting.R")

mark <- floor(seq(1, 600, length.out = 100))
sam_grid <- seq(0, 600, length.out = 100)
ori_dat <- dat
dat <- dat[mark,]
J <- length(sam_grid)

ggplot() + 
  geom_point(aes(x = dat[,1], y = dat[,2]), size = 0.6) + 
  labs(x = "X(cm)", y = "Z(cm)",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))
ggsave(paste0("Plot/data_1.pdf"), width = 4, height = 3, dpi = 300)

ggplot() + 
  geom_point(aes(x = seq(0, 6, length.out = nrow(dat)), y = dat[,1]), size = 0.6) + 
  labs(x = "Time(second)", y = "X(cm)",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))
ggsave(paste0("Plot/data_2.pdf"), width = 4, height = 3, dpi = 300)

ggplot() + 
  geom_point(aes(x = seq(0, 6, length.out = nrow(dat)), y = dat[,2]), size = 0.6) + 
  labs(x = "Time(second)", y = "Z(cm)",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))
ggsave(paste0("Plot/data_3.pdf"), width = 4, height = 3, dpi = 300)

# Parameter estimation
para_es <- para_est(dat, ord, para_length, int_upp, omega, para, sam_grid, diff = T, linear = T)

## Generalized smoothing approach
Gene_smo_res <- prof(int_upp, J, sam_grid, dat, para_es$Result_3)
  
# Plot
## Reconstruction
smo_dat <- para_es$smo_dat
init <- c(smo_dat[1,,1], smo_dat[1,,2])
time_grid <- seq(0, 600, length.out = 600)

gra_mat_res <- ode(y = init, times = time_grid, func = dynamic_equ, 
                   parms = para_es$Result_1, dim = p, omega = omega)[,2:3]
int_mat_res <- ode(y = init, times = time_grid, func = dynamic_equ, 
                   parms = para_es$Result_2, dim = p, omega = omega)[,2:3]
gre_mat_res <- ode(y = init, times = time_grid, func = dynamic_equ, 
                   parms = para_es$Result_3, dim = p, omega = omega)[,2:3]
pro_mat_res <- ode(y = init, times = time_grid, func = dynamic_equ, 
                   parms = Gene_smo_res$pars, dim = p, omega = omega)[,2:3]

dat_plot <- rbind(gra_mat_res, int_mat_res, gre_mat_res, pro_mat_res)
dat_plot[dat_plot < - 0.2] <- NA
dat_plot[dat_plot > 0.2] <- NA
dat_plot <- data.frame(dat_plot, c(rep("    Order 2 gradient matching", 600), rep("   Order 1 gradient matching", 600), 
                                   rep("  Green's matching", 600), rep(" GSA", 600)))
colnames(dat_plot) <- c("X", "Z", "Method")
dat_plot <- data.frame(dat_plot)

dat_p <- data.frame(rbind(ori_dat, 
                          ori_dat,
                          ori_dat,
                          ori_dat),
                    c(rep("    Order 2 gradient matching", 600), rep("   Order 1 gradient matching", 600), 
                      rep("  Green's matching", 600), rep(" GSA", 600)))
colnames(dat_p) <- c("X", "Z", "Method")

dat_obe <- ori_dat
dat_obe[-mark,] <- NA
dat_obe <- data.frame(rbind(dat_obe, 
                            dat_obe,
                            dat_obe,
                            dat_obe),
                    c(rep("    Order 2 gradient matching", 600), rep("   Order 1 gradient matching", 600), 
                      rep("  Green's matching", 600), rep(" GSA", 600)))
colnames(dat_obe) <- c("X", "Z", "Method")

dat_plot <- rbind(dat_plot, dat_p, dat_obe)
dat_plot <- data.frame(dat_plot, c(rep("Reconstruction", 2400), rep("Test data", 2400), rep("Fitting data", 2400)))
colnames(dat_plot) <- c("X", "Z", "Method", "L")
dat_plot$L <- factor(dat_plot$L, levels = c("Test data", "Reconstruction", "Fitting data"))

library(ggthemes)
ggplot(dat_plot) + 
  geom_point(aes(x = X, y = Z, color = L)) +
  facet_wrap(~Method, nrow = 2) + 
  theme_bw(base_family = "Times") +
  theme_tufte() +
  labs(x = "X(cm)", y = "Z(cm)",
       title = "",
       colour = "", fill = "", linetype = "") +
  scale_size_manual(values = c(1, 2, 2)) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        text = element_text(size = 17),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) 
  
ggsave(paste0("Plot/Fitting.pdf"), width = 9, height = 6.5, dpi = 300)

