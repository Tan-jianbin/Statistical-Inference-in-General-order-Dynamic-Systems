# Setting address
# setwd("~") # Set your address

# Loading model
source("Equation_discovery/Model.R")

# Running
## setting
J <- 50
rat <- 0.05

## Replication result
Result <- lapply(1:100, function(seed){
  sim_function(seed, p, J, para, rat, ord, int_upp, omega)
})

### True
dat_plot <- sapply(seq(-2, 2, length.out = 200), function(a){
  sapply(seq(-2, 2, length.out = 200), function(b){
    c <- con_state(matrix(c(a, b), nrow = 1), 0) %*% con_coef(para)
    d <- sqrt(b ^ 2 + c ^ 2)
    return(c(a, b, c, d))
  }, simplify = "array")
}, simplify = "array")

dat_plot <- apply(dat_plot, 1, function(a) a)
colnames(dat_plot) <- c("X", "DX", "D2X", "Length")
dat_plot_ture <- data.frame(dat_plot)

### Gradient matching
dat_plot <- sapply(seq(-2, 2, length.out = 200), function(a){
  sapply(seq(-2, 2, length.out = 200), function(b){
    c_sim <- sapply(1:length(Result), function(k){
      as.numeric(con_state(matrix(c(a, b), nrow = 1), 0) %*% Result[[k]]$para$Result_1$beta)
    })
    c <- mean(c_sim)
    d_sim <- sqrt(b ^ 2 + c_sim ^ 2)
    d <- mean(d_sim)
    return(c(a, b, c, d))
  }, simplify = "array")
}, simplify = "array")

dat_plot <- apply(dat_plot, 1, function(a) a)
colnames(dat_plot) <- c("X", "DX", "D2X", "Length")
dat_plot_grad <- data.frame(dat_plot)

### Green's matching
dat_plot <- sapply(seq(-2, 2, length.out = 200), function(a){
  sapply(seq(-2, 2, length.out = 200), function(b){
    c_sim <- sapply(1:length(Result), function(k){
      as.numeric(con_state(matrix(c(a, b), nrow = 1), 0) %*% (Result[[k]]$para$Result_2$beta)[1:7])
    })
    c <- mean(c_sim)
    d_sim <- sqrt(b ^ 2 + c_sim ^ 2)
    d <- mean(d_sim)
    return(c(a, b, c, d))
  }, simplify = "array")
}, simplify = "array")

dat_plot <- apply(dat_plot, 1, function(a) a)
colnames(dat_plot) <- c("X", "DX", "D2X", "Length")
dat_plot_gree <- data.frame(dat_plot)

dat_plot <- rbind(dat_plot_ture, dat_plot_grad, dat_plot_gree)
dat_plot <- data.frame(dat_plot, Type = c(rep(" True", 40000), rep("  Gradient matching", 40000), rep("Green's matching", 40000)))
dat_plot$Length[dat_plot$Length > 2.5] <- NA 

ggplot(dat_plot, aes(x = X, y = DX)) +
  geom_point(aes(color = Length)) +
  facet_wrap(~Type) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_colour_gradientn(
    limits = c(0, 2.5),
    # breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
    colours = c("red", 'orange', 'yellow','purple', "blue"),
    space = "Lab",
    na.value = "blue",
    guide = "colourbar",
    aesthetics = "colour"
  ) +
  labs(x = "X", y = "DX",
       title = "",
       colour = "", fill = "", linetype = "") +
  ggthemes::theme_tufte() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        text = element_text(size = 15),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))

ggsave(paste0("Plot/vector_field.pdf"), width = 10, height = 4, dpi = 300)
