########################################################################

# FUNCTIONS FOR THE MODEL

########################################################################

# Library
#---------------

library(doParallel)
library(tidyverse)
library(plotly)
library(patchwork)
library(RANN)

# Functions
#---------------

forest_new <- function(coefficients, density, n_reprod = 1) {
    N_layers = nrow(density)
    N_sp = ncol(density)

    # Extract coefficients for the given species
    coef <- coefficients[colnames(density), ]
    growth1 <- coef$growth1
    growth2 <- coef$growth2
    birth <- coef$birth
    mortality <- coef$mortality
    LC_growth <- coef$LC_growth
    LC_birth <- coef$LC_birth
    LC_mortality <- coef$LC_mortality

    # Calculate the modulator (i.e. list of the total cover density above (and containing) the layer of interest)
    Modulateur <- c()
    for(layer in 1:N_layers){
        Modulateur[layer] <- basal_area[layer] * sum(as.numeric(density[layer,]))
    }
    Modulateur <- Modulateur %>% rev()
    Modulateur <- cumsum(Modulateur) %>%  rev()

    # Create a new matrix for the new densities
    density_new <- matrix(0, nrow = N_layers, ncol = N_sp)
    for (sp in 1:N_sp) {
        density_new[1,sp] <- max(0,
            (density[1,sp] + birth[sp] * (1 - LC_birth[sp] * Modulateur[1]) +
            - density[1,sp] * mortality[sp] * (LC_mortality[sp] * Modulateur[2] + 1) +
            - (growth1[sp]/22.5) * density[1,sp] * (1 - LC_growth[sp] * Modulateur[2])) %>% as.numeric())
        
        density_new[2,sp] <- max(0,
            (density[2,sp] + 
            + (growth1[sp]/22.5) * density[1,sp] * (1 - LC_growth[sp] * Modulateur[2]) +
            - (growth2[sp]/45) * density[2,sp] * (1 - LC_growth[sp] * Modulateur[3]) +
            - density[2,sp] * mortality[sp] * (LC_mortality[sp] * Modulateur[3] + 1)) %>%
            as.numeric())
        
        density_new[3,sp] <- max(0,
            (density[3,sp] +
            + (growth2[sp]/45) * density[2,sp] * (1 - LC_growth[sp] * Modulateur[3]) +
            - density[3,sp] * mortality[sp]) %>%
            as.numeric())
    }
    colnames(density_new) <- colnames(density)
    rownames(density_new) <- rownames(density)
    return(density_new)
}

simul_forest <- function(coefficients, density, T = 200, control = NULL){
    # Error if you did not give enough basal area data for the number of layer
    if(length(basal_area) < nrow(density)){stop("basal_area must have the same length as the number of layer")}
    # Warning if you give too much basal area data for the number of layer
    if(length(basal_area) > nrow(density)){warning("basal_area is longer than the number of layer")}

    # Get the initial state of the forest and associated metrics
    forest <- cbind(density, layer = 1:nrow(density), time = 1, control = FALSE, data.frame(rbind(forest_metric(density))), extraction = NA)

    # Simulate the forest for T time step
    for (t in 2:T){
        if(t%%5 == 2 & !is.null(control)){ # Every 5 time step, apply the control beginning with the first time step
            old_Biomasse = data.frame(rbind(forest_metric(density)))[,2]
            density<-apply_control(density, control[[(t+3)/5]])
            extraction = old_Biomasse - data.frame(rbind(forest_metric(density)))[,2] # Calculate the extraction
            forest <- rbind(forest, cbind(density, layer = 1:nrow(density), time = t - 1, control = TRUE, data.frame(rbind(forest_metric(density))), extraction = extraction))
        }
        density <- forest_new(coefficients, density, n_reprod)
        forest <- rbind(forest, cbind(density, layer = 1:nrow(density), time = t, control = FALSE, data.frame(rbind(forest_metric(density))), extraction = NA))
    }
    return(forest)
}

simul_forest_fast <- function(param, density, T = 20){
    forest <- cbind(density, layer = 1:nrow(density), time = 1)
    for (t in 2:T){
        density <- forest_new(param, density, n_reprod = 1)
        forest <- rbind(forest, cbind(density, layer = 1:nrow(density), time = t))
    }
    return(data.frame(forest) %>% pivot_longer(cols = colnames(density), names_to = "species", values_to = "nb_trees"))
}

#' Forest_new_5
#' 
#' This function calculate the new density of a forest for 5 time step from the previous density
#' 
#' @param density A matrix of the density of each species in each layer
#' @param n_reprod The number of layer where reproduction is possible
#' @return A matrix of the new density of each species in each layer

forest_new_5 <- function(param, density, n_reprod = 1){
    Year = forest_new(param, density, n_reprod)
    for(i in 1:4){
        Year = forest_new(param, Year, n_reprod)
    }
    return(Year)
}

#' Apply_control
#' 
#' This function apply a control to a forest (control is defined as the number of trees left in each layer for each species)
#' 
#' @param density A matrix of the density of each species in each layer
#' @param control A matrix of the control to apply (NA if the layer is not controlled)
#' @return The forest density with the control applied

apply_control <- function(density, control){
    for(layer in 1:nrow(density)){
        for(sp in 1:ncol(density)){
            if(is.na(control[layer, sp]) | control[layer, sp] > density[layer, sp]){next}
            # if the control (nb of tree to leave) is possible then apply it
            density[layer, sp] <- control[layer, sp] %>% as.numeric()
        }
    }
    return(density)
}

#' Forest_metric
#' 
#' This function calculate the metric of a forest : number of trees, biomass, Shannon index
#' 
#' @param density A matrix of the density of each species in each layer
#' @return A data frame with the number of trees, biomass, Shannon index

forest_metric <- function(density){
    N_layers = nrow(density)
    N_tree = sum(density)
    Basal_area = sum(density * basal_area[1:N_layers])
    Vol_wood = sum(Volume[1:N_layers,colnames(density)] * density)
    Shannon_vert = -sum(rowSums(density)/N_tree * log2(rowSums(density)/N_tree))
    Shannon= -sum(colSums(density)/N_tree * log2(colSums(density)/N_tree))
    #ShannoN_layers = apply(density, 1, function(x) -sum(x/N_tree * log2(x/N_tree)))
    metric <- c(N_tree, Basal_area, Vol_wood, Shannon_vert, Shannon)
    return(metric)
}

#' Reform_forest
#' 
#' This function reform the forest data frame to be used in ggplot by pivot longer
#' 
#' @param forest A data frame with the density of each species in each layer for each time step
#' @param species_selection The species present in the forest
#' @return A data frame with the density of each species in each layer for each time step and the layer as a factor, (columns : time, layer, Species, Density)

reform_forest <- function(forest, species_selection){
    return(forest %>% data.frame() %>%
        pivot_longer(cols = species_selection, names_to = "species", values_to = "density") %>%
        mutate(layer = as.character(layer), layer = factor(layer, level = rev(unique(forest$layer)))))
}

#' Graph_forest
#' 
#' This function plot the forest data frame for multiple species and layer with ggplot

multi_sim <- function(min = 0, T = 200, N_layers = 2){
    # Simulate a lot of forest with different mixture and initial states
    all_sp <- list(c("Sapin", "Epicea"), c("Sapin", "Pin"), c("Sapin", "Bouleau"), c("Sapin", "Hetre"), c("Sapin", "Chene"), c("Epicea", "Pin"), c("Epicea", "Bouleau"), c("Epicea", "Hetre"),
    c("Epicea", "Chene"), c("Pin", "Bouleau"), c("Pin", "Hetre"), c("Pin", "Chene"), c("Bouleau", "Hetre"), c("Bouleau", "Chene"), c("Hetre", "Chene"))
    EI = matrix(min, nrow = N_layers, ncol = 2) %>% data.frame()
    
    forest <- data.frame()
    for (i in 1:length(all_sp)){
        init = EI
        colnames(init) <- all_sp[[i]]
        forest <- forest %>% 
            rbind(cbind(simul_forest(coefficients, init, T, control = NULL) %>% reform_forest(all_sp[[i]]), association1 = all_sp[[i]][1], association2 = all_sp[[i]][2]))
    }
    return(forest)
}

#' Simple_simul
#' 
#' This function simulate a forest for a given number of time step with a uniform initial state, it is used to make Simul_forest easier to use
#' 
#' @param unif_is The uniform initial state
#' @param N_layers The number of layer
#' @param species_selection The species present in the forest
#' @param n_reprod The number of layer where reproduction is possible
#' @param control A list of matrix of the density of each species in each layer for each time step where control is applied
#' @param T The number of time step
#' @return The forest simulated with the density of each species in each layer for each time step and the metric of the forest

Simple_simul <- function(param, unif_is, N_layers, species_selection, n_reprod = 1, control = NULL, T = 200){
    density<- matrix(rep(unif_is, N_layers * length(species_selection)), nrow = N_layers, ncol = length(species_selection)) %>% data.frame()
    colnames(density) <- species_selection
    forest <- simul_forest(param, density, T, control = control) %>% reform_forest(species_selection)
    forest <- forest %>% select(time, layer, species, density)
    return(forest)
}

#' Graph_forest
#' 
#' This function plot the forest data frame for multiple species and layer with ggplot
#' 
#' @param forest A data frame with the density of each species in each layer for each time step
#' @return A ggplot object

graph_forest <- function(forest){
    ggplot(forest, aes(x = time, y = density, color = species, linetype = factor(layer))) +
        geom_line() +
        theme_bw() +
        labs(x = "Time", y = "Density")
}

################
#---------------
# Necessary data
#---------------
################

#coefficients <- data.frame(
#    Species = c("Sapin", "Epicea", "Pin", "Bouleau", "Hetre", "Chene"),
#    growth = c(0.025, 0.018, 0.02, 0.025, 0.08),
#    birth = c(0.75, 0.75, 0.5, 0.75, 1.5),
#    mortality = c(0.0067, 0.005, 0.002, 0.006, 0.01),
#    LC_growth = c(0.0167, 0.0167, 0.02, 0.022, 0.03),
#    LC_birth = c(0.0125, 0.0125, 0.016, 0.014, 0.03),
#    LC_mortality = c(0.0008, 0.0008, 0.001, 0.001, 0.012)
#)
#rownames(coefficients) <- c("Sapin", "Epicea", "Pin", "Bouleau", "Hetre", "Chene")
#
#basal_area = c(0.013, 0.16, 0.18)
#
#L_sp <- c(3,3,3,3,2)
#
#Volume = data.frame(
#    Resineux = c(0.013, 0.16, 0.18), 
#    Hetre = c(0.013, 0.16, 0.18), 
#    Chene = c(0.013, 0.16, 0.18), 
#    Pin = c(0.013, 0.16, 0.18), 
#    Bouleau = c(0.013, 0.16, 0.18))
