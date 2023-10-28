source("function.R")

###################################
# data

coefficients <- data.frame(
    Species = c("Resineux", "Hetre", "Chene", "Pin", "Bouleau"),
    growth = c(0.025, 0.018, 0.02, 0.025, 0.03),
    birth = c(0.75, 0.75, 0.5, 0.75, 1.5),
    mortality = c(0.0067, 0.005, 0.003, 0.0067, 0.005),
    LC_growth = c(0.0167, 0.0167, 0.025, 0.02, 0.035),
    LC_birth = c(0.0125, 0.0125, 0.016, 0.016, 0.03),
    LC_mortality = c(0.0008, 0.0008, 0.001, 0.001, 0.002)
)


coefficients <- data.frame(
    Species = c("Resineux", "Hetre", "Chene", "Pin", "Bouleau"),
    growth = c(0.025, 0.018, 0.02, 0.025, 0.03),
    birth = c(0.75, 0.75, 0.5, 0.75, 1.5),
    mortality = c(0.0067, 0.005, 0.003, 0.0067, 0.005),
    LC_growth = c(0.0167, 0.0167, 0.025, 0.025, 0.03),
    LC_birth = c(0.0125, 0.0125, 0.016, 0.016, 0.03),
    LC_mortality = c(0.0008, 0.0008, 0.001, 0.001, 0.002))
rownames(coefficients) <- c("Resineux", "Hetre", "Chene", "Pin", "Bouleau")

multi_sim()
basal_area = c(0.013, 0.16)

Volume = data.frame(
    Resineux = c(0.013, 0.16, 0.18), 
    Hetre = c(0.013, 0.16, 0.18), 
    Chene = c(0.013, 0.16, 0.18), 
    Pin = c(0.013, 0.16, 0.18), 
    Bouleau = c(0.013, 0.13, 0.18))

###################################
# control exemple

c_Res = 20
c_Bou = 40
control = list(c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou),
    c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou),
    c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou),
    c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou),
    c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou), c(c_Res,c_Bou))

###################################
# Forest simulation exemple

forest1 <- forest_simul_list(100, c(70,90),c(100,90),control, c("Resineux","Bouleau"))

ggplot(forest1) +
    geom_line(aes(x = Time, y = Density, color = Species, linetype = Compartment)) +
    geom_line(aes(x = Time, y = Biomass, color = "Biomasse"), linetype = "dashed") +
    geom_line(aes(x = Time, y = Tot_prod, color = "Production cumulée"), linetype = "dashed") +
    geom_line(aes(x = Time, y = Shannon, color = "Diversité"), linetype = "dashed") +
    theme(legend.position = "none")+
    labs(x = "Temps", y = "Densité") +
    theme_bw() +
    ylim(0,300)

multi_sim(EI_min = c(20,20), EI_max = c(100,100))

###################################
# Viability kernel exemple

f = 50
x1 = seq(0, 300, f)
x2 = seq(0, 300, f)
fc = 50
u = seq(0,100,fc)

table <- forest_transition_table(c("Resineux","Bouleau"))
saveRDS(table, "table.rds")
readRDS("table.rds") -> table

table_v <- table_viability(table, sh = 0.6, biom = 10, prod = 20)
table_v %>% filter(Constraints == 1) %>% nrow()
table_v %>% filter(Viable == 1) %>% nrow()
table_v %>% filter(Viable_nn == 1) %>% nrow()

Shannon = seq(0,1,0.1)
Biomass = seq(0,200,10)
Extraction = seq(0,200,10)
sensibilite <- data.frame()
# tester la sensibilité aux différentes contraintes
for(s in Shannon){
    for(b in Biomass){
        for(p in Extraction){
            table_v <- table_viability(table, sh = s, biom = b, prod = p)
            sensibilite <- rbind(sensibilite,
                c(s,b,p,table_v %>% filter(Constraints == 1) %>% nrow(),
                table_v %>% filter(Viable == 1) %>% nrow()))
        }
    }
}

colnames(sensibilite) <- c("Shannon", "Biomass", "Extraction", "Constraints", "Viable")
short_sensibilite <- sensibilite %>%
    group_by(Shannon, Biomass, Extraction) %>%
    summarize(Constraints = sum(Constraints)/n(), Viable = sum(Viable)/n())
library(plotly)
plot_ly(short_sensibilite %>% filter(Constraints > 0), x = ~Shannon, y = ~Biomass, z = ~Extraction, color = ~Viable>0)
plot_ly(short_sensibilite %>% filter(Viable >0), x = ~Shannon, y = ~Biomass, z = ~Extraction, color = ~Viable) %>%
    layout(
    xaxis = list(range=c(0,1)), yaxis = list(range=c(0,200)), zaxis = list(range=c(0,200)))

head(short_sensibilite)

# Ajouter la colonne du nombre de possibilités viable à partir d'un control
table_v %>% group_by(x1, x2) %>% mutate(flexi = sum(Viable)) -> table_v

# séparer les colonens de double x1 et x2, ainsi que u en x1_1, x1_2 puis x2_1, x2_2 et u_1, u_2
table_v %>% separate(x1, into = c("x1_1", "x1_2"), sep = ",") %>%
    separate(x2, into = c("x2_1", "x2_2"), sep = ",") %>%
    separate(u, into = c("u_1", "u_2"), sep = ",") -> table1
table1 %>% separate(x1_nn, into = c("x1_1_nn", "x1_2_nn"), sep = ",") %>%
    separate(x2_nn, into = c("x2_1_nn", "x2_2_nn"), sep = ",") -> table1

table1$x1_1 <- gsub("[c\\(\\)]", "", table1$x1_1)
table1$x1_2 <- gsub("[c\\(\\)]", "", table1$x1_2)
table1$x2_1 <- gsub("[c\\(\\)]", "", table1$x2_1)
table1$x2_2 <- gsub("[c\\(\\)]", "", table1$x2_2)
table1$u_1 <- gsub("[c\\(\\)]", "", table1$u_1)
table1$u_2 <- gsub("[c\\(\\)]", "", table1$u_2)
table1$x1_1_nn <- gsub("[c\\(\\)]", "", table1$x1_1_nn)
table1$x1_2_nn <- gsub("[c\\(\\)]", "", table1$x1_2_nn)
table1$x2_1_nn <- gsub("[c\\(\\)]", "", table1$x2_1_nn)
table1$x2_2_nn <- gsub("[c\\(\\)]", "", table1$x2_2_nn)

table1 <- table1 %>% mutate(
    x1_1 = as.numeric(x1_1),
    x1_2 = as.numeric(x1_2),
    x2_1 = as.numeric(x2_1),
    x2_2 = as.numeric(x2_2),
    u_1 = as.numeric(u_1),
    u_2 = as.numeric(u_2),
    x1_1_nn = as.numeric(x1_1_nn),
    x1_2_nn = as.numeric(x1_2_nn),
    x2_1_nn = as.numeric(x2_1_nn),
    x2_2_nn = as.numeric(x2_2_nn)
)

# Plot une dimension de Viability kernel
ggplot(table1 %>% filter(x1_1 == 250, x1_2 == 100, u_1 == 0, u_2 == 0)) +
    geom_point(aes(x = x2_1, y = x2_2, color = Viable))

ggplot(table1 %>% filter(x1_1 == 250, x1_2 == 100) %>% group_by(x2_2, x2_1) %>% summarize(Viable = sum(Viable))) +
    geom_point(aes(x = x2_1, y = x2_2, color = Viable))

#Faire un itinéraire au hasard
traj <- table1 %>% filter(Viable == 1) %>% sample_n(1)
for(i in 1:50){
    traj <- traj %>% rbind(table1 %>% filter(x1_1 == traj$x1_1_nn[i], x2_1 == traj$x2_1_nn[i], x1_2 == traj$x1_2_nn[i], x2_2 == traj$x2_2_nn[i]) %>%
    sample_n(1))
}
traj

# Faire un itinéraire en maximisant flexi
traj <- table1 %>% filter(Viable == 1) %>% sample_n(1)
for(i in 1:100){
    traj <- traj %>% rbind(table1 %>% filter(x1_1 == traj$x1_1_nn[i], x2_1 == traj$x2_1_nn[i], x1_2 == traj$x1_2_nn[i], x2_2 == traj$x2_2_nn[i]) %>%
    filter(flexi == max(flexi)) %>% sample_n(1))
}
traj

control = traj %>% rowwise() %>% mutate(u = list(c(u_1, u_2))) %>% pull(u)

forest1 <- forest_simul_list(50, c(traj$x1_1[1],traj$x1_2[1]),c(traj$x2_1[1],traj$x2_2[1]),control, c("Resineux","Bouleau"))
ggplot(forest1) +
    geom_line(aes(x = Time, y = Density, color = Species, linetype = Compartment)) +
    #geom_line(aes(x = Time, y = Biomass, color = "Biomasse"), linetype = "dashed") +
    #geom_line(aes(x = Time, y = Tot_prod, color = "Production cumulée"), linetype = "dashed") +
    #geom_line(aes(x = Time, y = Shannon, color = "Diversité"), linetype = "dashed") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    theme(legend.position = "none")+
    labs(x = "Temps", y = "Densité") +
    theme_bw() +
    ylim(0,300)

ggplot(table1 %>% group_by(u_1, u_2) %>% summarize(Viable = sum(Viable)) %>% filter(Viable != 0)) +
    geom_point(aes(x = u_1, y = u_2, color = Viable))

table1 <-  table1 %>% 
    mutate(
        U_div = u_1/(u_1 + u_2),
        U_biomass = u_1 + u_2, # Pas sure sure pour ce calcul
        U_ext = (x2_1 - u_1) + (x2_2 - u_2)
)
ggplot(table1) +
    geom_bar(aes(x = U_ext, color = factor(Viable), fill = factor(U_div)), position = "fill")

ggplot(table1 %>% group_by(U_div, U_biomass) %>% summarize(freq_v = sum(Viable)/n())) +
    geom_point(aes(x = U_div, y = U_biomass, color = freq_v))

ggplot(table1 %>% group_by(U_div, U_ext) %>% summarize(freq_v = sum(Viable)/n())) +
    geom_point(aes(x = U_div, y = U_ext, color = freq_v))
