source("Fonctions.R")

forest <- Simple_simul(unif_is = 100, N_comp = 3, species_selection = c("Resineux", "Bouleau"))
graph_forest <- function(forest){
    ggplot(forest, aes(x = time, y = Density, color = Species, linetype = factor(Compartment))) +
        geom_line() +
        theme_bw() +
        labs(x = "Temps", y = "Densité") +
        ylim(0,400)
}

# Que pour deux espèces mais donne tous les graphes
multi_sim(N_comp = 3)

basal_area = c(0.013, 0.16, 0.18)
# Autant d'espèce que l'on veut donne 1 graphe
graph_forest(Simple_simul(unif_is = 10, N_comp = 3, species_selection = c("Resineux"), T = 500))
graph_forest(Simple_simul(unif_is = 10, N_comp = 3, species_selection = c("Bouleau"))) +
graph_forest(Simple_simul(unif_is = 10, N_comp = 3, species_selection = c("Hetre"))) +
graph_forest(Simple_simul(unif_is = 10, N_comp = 3, species_selection = c("Chene"))) +
graph_forest(Simple_simul(unif_is = 10, N_comp = 3, species_selection = c("Pin")))



density

density = matrix(c(20,20,10), nrow = 3, ncol = 1)
colnames(density) = c("Resineux")
N_comp = 3
density <- forest_new(density)
density


N_comp = nrow(density)
N_sp = ncol(density)

# Extract coefficients for the given species
coef <- coefficients[colnames(density), ]
growth <- coef$growth
birth <- coef$birth
mortality <- coef$mortality
LC_growth <- coef$LC_growth
LC_birth <- coef$LC_birth
LC_mortality <- coef$LC_mortality

# Calculate the modulator
Modulateur <- c()
for(comp in 1:N_comp){
    Modulateur[comp] <- basal_area[comp] * sum(as.numeric(density[comp,]))
}
Modulateur <- Modulateur %>% rev()
Modulateur <- cumsum(Modulateur) %>%  rev()
# Create a new matrix for the new densities
density_new <- matrix(0, nrow = N_comp, ncol = N_sp)
for (sp in 1:N_sp) {
    strate_reproduction <- sum(basal_area[(N_comp-n_reprod + 1):N_comp] * density[(N_comp-n_reprod +1):N_comp,sp])
    Birth = as.numeric(birth[sp] * strate_reproduction * (1 - LC_birth[sp] * Modulateur[1]))
    Growth = as.numeric(- density[1,sp] * (LC_mortality[sp] * Modulateur[2] + mortality[sp]))
    Mortality = as.numeric(- growth[sp] * density[1,sp] * (1 - LC_growth[sp] * Modulateur[2]))
    density_new[1,sp] <- max(0,
        (density[1,sp] + 
        + birth[sp] * strate_reproduction * (1 - LC_birth[sp] * Modulateur[1]) +
        - density[1,sp] * (LC_mortality[sp] * Modulateur[2] + mortality[sp]) +
        - growth[sp] * density[1,sp] * (1 - LC_growth[sp] * Modulateur[2])) %>% as.numeric())
    print(c(Birth, Growth, Mortality, strate_reproduction))
    if(N_comp > 2){
    for(comp in 2:(N_comp - 1)){
        density_new[comp,sp] <- max(0,
            (density[comp,sp] + 
            + growth[sp] * density[comp-1,sp] * (1 - LC_growth[sp] * Modulateur[comp]) +
            - growth[sp] * density[comp,sp] * (1 - LC_growth[sp] * Modulateur[comp + 1]) +
            - density[comp,sp] * (LC_mortality[sp] * Modulateur[comp + 1] + mortality[sp])) %>%
            as.numeric())
    }}
    density_new[N_comp,sp] <- max(0,
        (density[N_comp,sp] + 
        + growth[sp] * density[N_comp-1,sp] * (1 - LC_growth[sp] * Modulateur[N_comp]) +
        - density[N_comp,sp] * mortality[sp]) %>%
        as.numeric())
}
colnames(density_new) <- colnames(density)
rownames(density_new) <- rownames(density)