########################################################################

# FUNCTIONS FOR THE MODEL AND VIABILITY

########################################################################

################
#---------------
# Library
#---------------
################

library(doParallel)
library(tidyverse)
library(plotly)
library(patchwork)
library(RANN)

################
#---------------
# Necessary data
#---------------
################

coefficients <- data.frame(
    Species = c("Resineux", "Hetre", "Chene", "Pin", "Bouleau"),
    growth = c(0.025, 0.018, 0.02, 0.025, 0.08),
    birth = c(0.75, 0.75, 0.5, 0.75, 1.5),
    mortality = c(0.0067, 0.005, 0.002, 0.006, 0.01),
    LC_growth = c(0.0167, 0.0167, 0.02, 0.022, 0.03),
    LC_birth = c(0.0125, 0.0125, 0.016, 0.014, 0.03),
    LC_mortality = c(0.0008, 0.0008, 0.001, 0.001, 0.012)
)
rownames(coefficients) <- c("Resineux", "Hetre", "Chene", "Pin", "Bouleau")

basal_area = c(0.013, 0.16, 0.18)

Volume = data.frame(
    Resineux = c(0.013, 0.16, 0.18), 
    Hetre = c(0.013, 0.16, 0.18), 
    Chene = c(0.013, 0.16, 0.18), 
    Pin = c(0.013, 0.16, 0.18), 
    Bouleau = c(0.013, 0.16, 0.18))

################
#---------------
# Functions
#---------------
################

# Function for the model and simulation
#--------------------------------------
########################################

#' Forest_new
#' 
#' This function calculate the new density of a forest for one time step from the previous density
#' 
#' @param density A matrix of the density of each species in each layer
#' @param n_reprod The number of layer where reproduction is possible
#' @return A matrix of the new density of each species in each layer

forest_new <- function(density, n_reprod = 1) {
    N_layers = nrow(density)
    N_sp = ncol(density)

    # Extract coefficients for the given species
    coef <- coefficients[colnames(density), ]
    growth <- coef$growth
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

    # Calculate the new density for each species in each layer
    for (sp in 1:N_sp) {
        # Total density of the cover where reproduction is possible
        layer_reproduction <- sum(basal_area[(N_layers-n_reprod + 1):N_layers] * density[(N_layers-n_reprod +1):N_layers,sp])

        density_new[1,sp] <- max(0,
            (density[1,sp] + 
            # Birth
            + birth[sp] * layer_reproduction * (1 - LC_birth[sp] * Modulateur[1]) +
            # Mortality
            - density[1,sp] * (LC_mortality[sp] * Modulateur[2] + mortality[sp]) +
            # Growth from this layer to the next
            - growth[sp] * density[1,sp] * (1 - LC_growth[sp] * Modulateur[2])) %>% as.numeric())
        if(N_layers > 2){
        for(layer in 2:(N_layers - 1)){
            density_new[layer,sp] <- max(0,
                (density[layer,sp] + 
                # Growth from the previous layer
                + growth[sp] * density[layer-1,sp] * (1 - LC_growth[sp] * Modulateur[layer]) +
                # Growth to the next layer
                - growth[sp] * density[layer,sp] * (1 - LC_growth[sp] * Modulateur[layer + 1]) +
                # Mortality
                - density[layer,sp] * (LC_mortality[sp] * Modulateur[layer + 1] + mortality[sp])) %>%
                as.numeric())
        }}
        density_new[N_layers,sp] <- max(0,
            (density[N_layers,sp] + 
            # Growth from the previous layer
            + growth[sp] * density[N_layers-1,sp] * (1 - LC_growth[sp] * Modulateur[N_layers]) +
            # Mortality
            - density[N_layers,sp] * mortality[sp]) %>%
            as.numeric())
    }
    colnames(density_new) <- colnames(density)
    rownames(density_new) <- rownames(density)
    return(density_new)
}

#' Forest_new_5
#' 
#' This function calculate the new density of a forest for 5 time step from the previous density
#' 
#' @param density A matrix of the density of each species in each layer
#' @param n_reprod The number of layer where reproduction is possible
#' @return A matrix of the new density of each species in each layer

forest_new_5 <- function(density, n_reprod = 1){
    Year = forest_new(density, n_reprod)
    for(i in 1:4){
        Year = forest_new(Year, n_reprod)
    }
    return(Year)
}

#' Simul_forest
#' 
#' This function simulate a forest for a given number of time step
#' 
#' @param density A matrix of the initial density of each species in each layer
#' @param T The number of time step
#' @param control A list of matrix of the density of each species in each layer for each time step where control is applied
#' @return The forest simulated with the density of each species in each layer for each time step and the metric of the forest

simul_forest <- function(density, T = 100, control = NULL){
    # Error if you did not give enough basal area data for the number of layer
    if(length(basal_area) < nrow(density)){stop("basal_area must have the same length as the number of layer")}
    # Warning if you give too much basal area data for the number of layer
    if(length(basal_area) > nrow(density)){warning("basal_area is longer than the number of layer")}

    # Get the initial state of the forest and associated metrics
    forest <- cbind(density, layer = 1:nrow(density), time = 1, control = FALSE, data.frame(forest_metric(density)), Extraction = NA)

    # Simulate the forest for T time step
    for (t in 2:T){
        if(t%%5 == 2 & !is.null(control)){ # Every 5 time step, apply the control beginning with the first time step
            old_Biomasse = forest_metric(density)$Biomass
            density<-apply_control(density, control[[(t+3)/5]])
            extraction = old_Biomasse - forest_metric(density)$Biomass # Calculate the extraction
            forest <- rbind(forest, cbind(density, layer = 1:nrow(density), time = t - 1, control = TRUE, data.frame(forest_metric(density)), Extraction = extraction))
        }
        density <- forest_new(density, n_reprod = 1)
        forest <- rbind(forest, cbind(density, layer = 1:nrow(density), time = t, control = FALSE, data.frame(forest_metric(density)), Extraction = NA))
    }
    return(forest)
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
    Biomass = sum(Volume[1:N_layers,colnames(density)] * density)
    Shannon_vert = -sum(rowSums(density)/N_tree * log2(rowSums(density)/N_tree))
    Shannon= -sum(colSums(density)/N_tree * log2(colSums(density)/N_tree))
    ShannoN_layers = apply(density, 1, function(x) -sum(x/N_tree * log2(x/N_tree)))
    metric <- data.frame(N_tree, Biomass, Shannon_vert, Shannon)
    # change NA to 0
    metric[is.na(metric)] <- 0
    colnames(metric) <- c("N_tree", "Biomass", "Shannon_vert", "Shannon")
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
        pivot_longer(cols = species_selection, names_to = "Species", values_to = "Density") %>%
        mutate(layer = as.character(layer), layer = factor(layer, level = rev(unique(forest$layer)))))
}

#' Graph_forest
#' 
#' This function plot the forest data frame for multiple species and layer with ggplot

multi_sim <- function(min = 20, max = 100, T = 200, N_layers = 2){
    # Simulate a lot of forest with different mixture and initial states
    all_sp <- list(c("Resineux","Bouleau"), c("Resineux","Hetre"), c("Resineux","Chene"), c("Bouleau","Hetre"),
        c("Bouleau","Chene"), c("Bouleau","Pin"), c("Hetre","Chene"), c("Hetre","Pin"), c("Chene","Pin"))
    EI_min = matrix(min, nrow = N_layers, ncol = 2) %>% data.frame()
    EI_max = matrix(max, nrow = N_layers, ncol = 2) %>% data.frame()
    all_EI <- list(EI_min, EI_max)
    
    forest <- data.frame()
    for (i in 1:9){
        for(j in 1:2){
        init = all_EI[[j]]
        colnames(init) <- all_sp[[i]]
        forest <- forest %>% 
            rbind(cbind(simul_forest(init, T, control = NULL) %>% reform_forest(all_sp[[i]]), j, association1 = all_sp[[i]][1], association2 = all_sp[[i]][2]))
        }
    }
    ggplot(forest, aes(x = time, y = Density, color = Species, linetype = factor(layer))) +
        geom_line() +
        #scale_linetype_manual(values = c("solid", "dotted")) +
        theme_bw() +
        facet_grid(j~ paste0(association1, " - ", association2)) +
        labs(x = "Temps", y = "Densité") +
        ylim(0,300)
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

Simple_simul <- function(unif_is, N_layers, species_selection, n_reprod = 1, control = NULL, T = 200){
    density<- matrix(rep(unif_is, N_layers * length(species_selection)), nrow = N_layers, ncol = length(species_selection)) %>% data.frame()
    colnames(density) <- species_selection
    forest <- simul_forest(density, T, control = control) %>% reform_forest(species_selection)
    return(forest)
}

# Function for the viability
#--------------------------------------
########################################

transition_table <- function(N_layers, species_selection, N_cut = 1, n_reprod = 1, max_tree = 100, freq = 50){
    N_sp = length(species_selection)
    init <- seq(0, max_tree, freq)

    table <- expand.grid(rep(list(init), N_cut * N_sp + N_layers * N_sp))
    c = 0
    for(sp in 1:N_sp){
        for(cut in 1:N_cut){
            c = c + 1
            layer  = N_layers - N_cut + cut
            table$test <- table[,paste0("Var", c)] >=  table[,paste0("Var", N_cut * N_sp + (sp-1) * N_layers + layer)]
            table <- table %>% filter(test) %>% select(-test)
    }}

    # calculer les métriques
    metric <- metric_parallel(table[,columns_etat_initial(table, N_layers, species_selection, N_cut)], N_layers, species_selection)
    prod <- as.numeric(metric[,2]) - as.numeric(metric_parallel(table[,(N_cut * N_sp + 1):(N_cut* N_sp + N_sp * N_layers)], N_layers, species_selection)[,2]) 
    table <- table %>% cbind(fill_parallel(table[,-c(1: N_cut * N_layers)], N_layers, species_selection))
    metric <- cbind(metric, prod)
    colnames(metric) <- c("N_tree", "Biomass", "Shannon_vert", "Shannon", "ShannoN_layers", "Extraction")
    table <- table %>% cbind(metric)
    # plus proche voisin
    table$etat_initial <-  table[,columns_etat_initial(table, N_layers, species_selection, N_cut)] %>% group_indices(across(everything()), .add = TRUE)
    idx <- nn2(table[,columns_etat_initial(table, N_layers, species_selection, N_cut)], table[,(N_cut * N_sp + N_sp * N_layers + 1): (N_cut * N_sp + N_sp * N_layers + N_sp * N_layers)], k = 1)$nn.idx
    table$etat_final <- table[idx,"etat_initial"]
    return(table)
}

columns_etat_initial <- function(table, N_layers, species_selection, N_cut){
    N_sp = length(species_selection)
    columns <- c()
    for(sp in 1:N_sp){
        columns <- columns %>% append(c((N_cut * N_sp + 1 + (sp-1) * N_layers):(N_cut * N_sp + (sp-1) * N_layers + N_layers - N_cut)))
        #columns <- columns %>% append(c((N_cut * N_sp + (sp-1) * N_layers + N_layers - N_cut + 1): (N_cut * N_sp + (sp-1) * N_layers + N_layers - N_cut + 1 + N_cut - 1)))
        columns <- columns %>% append(c(((sp-1) * N_cut + 1): ((sp-1) * N_cut + N_cut)))
    }
    return(columns)
}

fill_parallel <- function(x, N_layers, species_selection){
    Ncpus <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(Ncpus)
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=1:nrow(x), .combine='rbind', .packages = c("tidyverse"),
      .export = c("forest_new_5", "forest_new", "forest_new_table", "coefficients", "basal_area")) %dopar%{
        return(forest_new_table(x[i,], N_layers, species_selection))
    }
  parallel::stopCluster(cl)
  return(res)
}

metric_parallel <- function(x, N_layers, species_selection){
    Ncpus <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(Ncpus)
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=1:nrow(x), .combine='rbind', .packages = c("tidyverse"),
      .export = c("forest_metric_table", "forest_metric", "coefficients", "basal_area", "Volume")) %dopar%{
        return(forest_metric_table(x[i,], N_layers, species_selection))
    }
  parallel::stopCluster(cl)
  colnames(res) <- c("N_tree", "Biomass", "Shannon_vert", "Shannon", "Shannon_layers")
  return(res)
}

forest_new_table <- function(row_table, N_layers, species_selection){
    N_sp = length(species_selection)
    density = matrix(row_table %>% as.numeric(), nrow = N_layers, ncol = N_sp) %>% data.frame()
    colnames(density) = species_selection
    density = forest_new_5(density, n_reprod = 1) %>% data.frame() %>% unlist() %>% as.numeric()
    return(density)
}

forest_metric_table <- function(row_table, N_layers, species_selection){
    N_sp = length(species_selection)
    density = matrix(row_table %>% as.numeric(), nrow = N_layers, ncol = N_sp) %>% data.frame()
    colnames(density) = species_selection
    metric = forest_metric(density)
    return(metric)
}

table_viability <- function(table, biomasse, shannon, ext, shannon_vert){
    # choix de contraintes
    table$Viable <- as.numeric(table$Biomass > biomasse & table$Shannon > shannon & table$Extraction > ext & table$Shannon_vert > shannon_vert)
    
    V <- table %>% select(etat_initial, etat_final, Viable)
    V$Constraints  = V$Viable
    for(i in 1:20){
        Viability <- V %>% select(etat_initial,Viable) %>% group_by(etat_initial) %>% summarize(Viable = sum(Viable)) %>%
            mutate(Viable = ifelse(Viable == 0, 0, 1))
        colnames(Viability) <- c("etat_final", "Viable_nn")
        V <- V %>% left_join(Viability, by = c("etat_final"))
        V <- V %>% mutate(Viable_new = Viable * Viable_nn)
        if(all(V$Viable == V$Viable_new)){break} 
        V <- V %>% select(-Viable_nn)
        V$Viable = V$Viable_new
    }
    print(paste("calculated whith : ", i, "iterations"))
    return(V %>% select(Constraints, Viable, Viable_nn))
}

columns_control <- function(table, N_layers, species_selection, N_cut, i){
    N_sp = length(species_selection)
    columns <- c()
    for(sp in 1:N_sp){
        columns <- columns %>% append(rep(NA, N_layers - N_cut))
        columns <- columns %>% append(table[i,((sp-1) * N_cut + 1): ((sp-1) * N_cut + N_cut)])
    }
    return(columns)
}