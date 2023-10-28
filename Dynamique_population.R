source("Fonctions.R")

plot_dynamic(unif_is = 20, N_comp = 2, species_selection = c("Resineux", "Bouleau"), n_reprod = 1)

multi_sim()

density <- matrix(c(20,20,90,20,70,10,20,20,80), nrow = 3, ncol = 3) %>% data.frame()
colnames(density) <- c("Resineux", "Hetre", "Chene")
control <- matrix(c(NA,NA,20,NA,NA,70,NA,NA,100), nrow = 3, ncol = 3) %>% data.frame()
colnames(control) <- c("Resineux", "Hetre", "Chene")
# liste de 10 fois les controles
control <- list(control, control, control, control, control, control, control, control, control, control)
forest <- simul_forest(density, T = 50, control = control)