library(tidyverse)
library(patchwork)

# parameters
g <- 0.025 # growth rate
m <- 0.017 # mortality rate
b <- 0.75 # birth rate
cg <- 0.0067 # competition growth
cb <- 0.0125 # competition birth
cm <- 0.0008 # competition mortality
g1 <- 0.013
g2 <- 0.16
v1 <- 0.066
v2 <- 2.29

x1_new <- function(x1,x2,u1){
    x = x1 +
    -g * x1 * (1 - cg * g2 * x2) + # growth
    b * g2 * x2 * (1 - cb * (g1 * x1 + g2 * x2)) + # birth
    - x1 * (cm * g2 * x2 + m) + # mortality
    - u1 * x1 # harvesting
    return(max(0,x))
}

x2_new <- function(x1,x2,u2){
    x = x2 +
    g * x1 * (1 - cg * g2 * x2) + # growth
    - m * x2  + # mortality
    - u2 * x2 # harvesting
    return(max(0,x))
}

forest_simul <- function(T,x1,x2,u = rep(1, 2*T/5)){
    forest <- data.frame(x1,x2, Prod = 0, u1 = 0, u2 = 0)
    Prod <- c(0)
    for (t in 1:T){# & u[(t+4)/5] * x1 * v1 + u[T/5 + (t+4)/5] * x2 * v2 > 30){
        h1 = 0; h2 = 0; p = 0
        if (t%%5 == 1){
            h1 = u[(t+4)/5]
            h2 = u[T/5 + (t+4)/5]
            p = h2 * x2 * v2 + h1 * x1 * v1}
        x1_new = x1_new(x1,x2,h1)
        x2 = x2_new(x1,x2,h2)
        x1 = x1_new
        forest <- forest %>% rbind(c(x1,x2,p, h1, h2))
    }
    # adding everything I want to study
    forest <- forest %>% 
        mutate(
            H = ifelse(x1*x2 == 0, 0, -((x1/(x1+x2)) * log2(x1/(x1+x2)) + (x2/(x1+x2)) * log2(x2/(x1+x2)))),
            t = 1:nrow(forest),
            Biomasse = x1 + x2)
    return(forest)
}

Optim_forest <- function(obj, T, x1, x2){

    forest_obj <- function(u){
    forest <- forest_simul(T = T, x1 = x1, x2 = x2, u)
    return(mean(forest[,obj]))
    }

    solution <- optim(rep(0,2 * T/5), forest_obj, lower = rep(0,2 * T/5),upper = rep(1,2 * T/5),method = "L-BFGS-B",  control = list(fnscale = -1))
    return(forest_simul(T = T, x1 = x1, x2 = x2, u = solution$par))
}

sol <- function(u1, u2){
    sol_x1 = c()
    sol_x2 = c()
    a_sol = - cg * g2^3 * b * cb * g
    b_sol = g2 * (m + u2) * (-g * cg + b * cb * g1 + cm) + g * b * g2^2 * (cg + cb)
    c_sol = (m + u2) * (g + m + u1) - g2 * b * g
    delta = b_sol^2 - 4 * a_sol * c_sol
    sol_x2[1] = (-b_sol + sqrt(delta)) / (2 * a_sol)
    sol_x2[2] = (-b_sol - sqrt(delta)) / (2 * a_sol)
    sol_x1[1] = (sol_x2[1] * (m + u2)) / (g - g * cg * g2 * sol_x2[1])
    sol_x1[2] = (sol_x2[2] * (m + u2)) / (g - g * cg * g2 * sol_x2[2])
    sol = c()
    for (i in 1:2){
        if(sol_x1[i]>=0 & sol_x2[i]>=0){
            sol <- append(sol, c(sol_x1[i], sol_x2[i]))
        }
    }
    if(length(sol) == 0){return(c(0,0))}
    return(sol)
}

f_H <- function(x1,x2){
    return(-((x1/(x1+x2)) * log2(x1/(x1+x2)) + (x2/(x1+x2)) * log2(x2/(x1+x2))))
}

f_Biomasse <- function(x1,x2){
    return(v1 * x1 + v2 * x2)
}

f_Prod <- function(x1,x2){
    return(0.1 * v1 * x1 + 0.05 * v2 * x2)
}

forest_5 <- function(x1,x2,u1,u2){
    Year = c(x1_new(x1,x2,u1), x2_new(x1,x2,u2))
    for(i in 1:4){
        Year = c(x1_new(Year[1], Year[2],0), x2_new(Year[1], Year[2],0))
    }
    return(Year)
}

##############################
# Model avec plusieurs espèces
##############################

# population décrit toujours avec 2 étages x1 et x2 mais avec N_sp espèces, on aura donc length(x1) = length(x2) = N_sp
# pour avoir l'état T + 1 de la foret il faudra calculer le modulateur de croissance, de mortalité et de naissace
# qui sera le même poour tout le monde : somme de x2_sp * g_sp = Modulateur (un réel)
# ensuite on pourra calculer x1_sp et x2_sp pour chaque espèce

#Tableau des données d'espèces :

forest_new_list <- function(x1,x2,u, sp){
    x1_new = c()
    x2_new = c()
    Mod_2 = sum(x2 * g2)
    Mod_1 = sum(x1 * g1)
    coef <- coefficients[sp,]
    g <- coef$g; b <- coef$b; m <- coef$m
    cg <- coef$cg; cb <- coef$cb; cm <- coef$cm
    for (i in 1:length(x1)){
        x1_new = append(x1_new, 
        max(0,x1[i] + -g[i] * x1[i] * (1 - cg * Mod_2) + # growth
        b[i] * g2[i] * x2[i] * (1 - cb[i] * (Mod_1 + Mod_2)) + # birth
        - x1[i] * (cm * Mod_2 + m[i]))) # mortality
        x2_new = append(x2_new, 
        max(0, x2[i] + g[i] * x1[i] * (1 - cg * Mod_2) + # growth
        - m[i] * x2[i]  + # mortality
        - u[i] * x2[i])) # harvesting
    }
    return(list(x1_new,x2_new))
}

forest_simul_list <- function(T,x1,x2,u, SP){
    forest <- data.frame(x1,x2,Prod = 0, u = 0, t = 1)
    forest$sp <- SP
    colnames(forest) <- c("x1","x2","Prod","u","t","sp")
    N_sp = length(x1)
    for (t in 2:T){
        h = c(0,0)
        #if (t%%5 == 1){
        #    h = u[[(t+4)/5]]}
        new_forest = forest_new_list(as.numeric(tail(forest$x1, N_sp)),as.numeric(tail(forest$x2, N_sp)),h, SP)
        new_forest = cbind(new_forest[[1]], new_forest[[2]], 0, 0 , t, SP)
        colnames(new_forest) <- c("x1","x2","Prod","u","t","sp")
        forest <- forest %>% rbind(new_forest)
    }
    # adding everything I want to study
    #forest <- forest %>% 
    #    mutate(
    #        H = ifelse(x1*x2 == 0, 0, -((x1/(x1+x2)) * log2(x1/(x1+x2)) + (x2/(x1+x2)) * log2(x2/(x1+x2)))),
    #        t = 1:nrow(forest),
    #        Biomasse = x1 + x2)
    # renommer les facteurs espèce avec le nom des espèces SP, remplacer 1 et 2 par les noms des espèces
    #forest <- forest %>% mutate(sp = factor(sp, levels = 1:length(sp), labels = sp))
        forest <- forest %>% mutate(
        x1 = as.numeric(x1), x2 = as.numeric(x2), Prod = as.numeric(Prod),
        u = as.numeric(u), t = as.numeric(t), sp = as.factor(sp)
    )
    forest <- forest %>% pivot_longer(cols = c(x1,x2), names_to = "compartiment", values_to = "density")
    return(forest)
}

all_sp <- list(c("Resineux","Bouleau"), c("Resineux","Hetre"), c("Resineux","Chene"), c("Bouleau","Hetre"),
    c("Bouleau","Chene"), c("Bouleau","Pin"), c("Hetre","Chene"), c("Hetre","Pin"), c("Chene","Pin"))
all_EI <- list(c(50,50), c(100, 100))

forest <- data.frame()
for (i in 1:9){
    for(j in 1:2){
    forest <- forest %>% 
        rbind(
            cbind(forest_simul_list(200, all_EI[[j]], all_EI[[j]],u,all_sp[[i]]), paste(all_sp[[i]][1],"-", all_sp[[i]][2]), j))
    }
}

colnames(forest) <- c("Prod","control","time","espece","compartiment","density", "association", "EI")

coefficients <- data.frame(
    g = c(0.025, 0.018, 0.02, 0.025, 0.03),
    b = c(0.75, 0.75, 0.5, 0.75, 1.5),
    m = c(0.0067, 0.005, 0.003, 0.0067, 0.005),
    cg = c(0.0167, 0.0167, 0.025, 0.025, 0.035),
    cb = c(0.0125, 0.0125, 0.016, 0.016, 0.03),
    cm = c(0.0008, 0.0008, 0.001, 0.001, 0.002)
)
rownames(coefficients) <- c("Resineux", "Hetre", "Chene", "Pin", "Bouleau")

forest1 <- forest_simul_list(100, c(200,200),c(200,200),u,c("Resineux","Bouleau"))

ggplot(forest) +
    geom_line(aes(x = time, y = density, color = espece, linetype = compartiment)) +
    theme(legend.position = "none")+
    labs(x = "Temps", y = "Densité") +
    theme_bw() +
    ylim(0,300) +
    facet_grid(EI ~ association)
