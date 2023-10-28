source("Fonctions.R")
library(plotly)

table <- readRDS("Data/table_long.RDS")
sensibilite <- readRDS("Data/sensibilite_long.RDS")

# tester le nombre d'états/contraintes viable par contraintes
short_sensibilite <- sensibilite %>%
    group_by(Shannon, Biomass, Extraction) %>%
    summarize(Constraints = sum(Constraints)/n(), Viable = sum(Viable)/n())
plot_ly(short_sensibilite, x = ~Shannon, y = ~Biomass, z = ~Extraction, color = ~Viable>0)
fig <- plot_ly(short_sensibilite, x = ~Shannon, y = ~Biomass, z = ~Extraction, color = ~Viable,
colorscale = 'Turbo') %>%
    layout(
    xaxis = list(range=c(0,1)), yaxis = list(range=c(0,200)), zaxis = list(range=c(0,200)))

# export plotly as html figures
htmlwidgets::saveWidget(
                widget = fig, #the plotly object
                file = "figure.html", #the path & file name
                selfcontained = TRUE #creates a single html file
                )

# Choisir une trajectoire au hasard
table_v <- cbind(table,table_viability(table, 10, 0.9, 3, 0.5)) %>% data.frame()
traj <- table_v %>% filter(Viable == 1) %>% sample_n(1)
for(i in 1:40){
    traj <- traj %>% rbind(table_v %>% filter(etat_initial == traj$etat_final & Viable == 1) %>%
    sample_n(1))
}
# Pour faire la vraie simul : a partir de traj
N_comp = 2
N_sp = 2
N_cut = 1
species_selection = c("Resineux","Bouleau")
density <- matrix(traj[1,columns_etat_initial(traj,N_comp, species_selection, N_cut)] %>% as.numeric(), nrow = N_comp, ncol = N_sp) %>% data.frame()
colnames(density) <- species_selection
control = list()
for(i in 1:nrow(traj)){
    control[[i]] <- matrix(columns_control(table,N_comp, species_selection, N_cut, i) %>% as.numeric(), nrow = N_comp, ncol = N_sp) %>% data.frame()
    colnames(control[[i]]) <- species_selection
}

forest <- simul_forest(density, T = 200, control = control) %>% reform_forest(species_selection)
ggplot(forest, aes(x = time, y = Density, color = Species, linetype = factor(Compartment))) +
    geom_line() +
    theme_bw() +
    labs(x = "Temps", y = "Densité") +
    geom_line(aes(y = Biomass))
View(head(forest))

traj$time = 0:(nrow(traj)-1) *5 +1

ggplot() +
    geom_line(data = traj, aes(x = as.numeric(time), y = as.numeric(Biomass)), color = "red") +
    geom_line(data = traj, aes(x = as.numeric(time), y = as.numeric(Extraction)), color = "blue") +
    geom_line(data = traj, aes(x = as.numeric(time), y = as.numeric(Shannon) * 10), color = "green") +
    geom_line(data = traj, aes(x = as.numeric(time), y = as.numeric(Shannon_vert) * 10), color = "yellow") +
    theme_bw() +
    labs(x = "Temps", y = "Biomasse")

table_v1 <- table_v %>% group_by(etat_initial) %>% mutate(Viability_state = sum(Viable),
    Viable = ifelse(Viable == 0, NA, Viable), Viability_state = ifelse(Viability_state == 0, NA, Viability_state))

plot_ly(table_v1 %>% filter(Var5 == 50), x = ~Var1, y = ~Var2, z = ~Var3, color = ~Viability_state, marker = list(size = 20), mode = 'markers') %>%
    layout(title = "Bouleau_1 = 50",scene = list(xaxis = list(title= "Resineux_2", range = c(0,100)), yaxis = list(title= "Bouleau_2", range = c(0,100)), zaxis = list(title= "Resineux_1", range = c(0,100))))