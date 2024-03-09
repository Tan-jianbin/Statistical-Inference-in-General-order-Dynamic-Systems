# Function
source("Function.R")

# Simulation
sfInit(parallel = T, cpus = 40)
sfExport("p", "J", "para", "ord", "para_length", "omega", "int_upp", "rat", "diff", "linear")
sfSource("Function.R")
sfExport("con_coef", "con_state", "dat_mean_gen", "force_func", "dynamic_equ", "jocobi")
Result <- sfLapply(seq(1, 200, length.out = 100), sim_function,
                   p = p, J = J,
                   rat = rat, para = para, ord = ord, linear = linear,
                   para_length = para_length, int_upp = int_upp,
                   omega = omega, diff = diff)
sfStop()
save(Result, file = paste0("Result/Result_", J, "_", rat, "_", mol_lab, ".rda"), version = 2)

