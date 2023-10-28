library(tidyverse)
library(patchwork)
library(RANN)

options(dplyr.summarise.inform = FALSE)

# Function

forest_new_list <- function(density_1, density_2, sp) {
    # forest_new_list function
    # Input:
    #   density_1: Initial density for compartment 1
    #   density_2: Initial density for compartment 2
    #   sp: Species index
    # Output:
    #   List containing new densities for compartment 1 and compartment 2

    # Initialize new density vectors
    density_1_new <- numeric()
    density_2_new <- numeric()
    Basal_area_1 <- 0.013
    Basal_area_2 <- 0.16

    # Calculate model parameters
    Mod_2 <- sum(density_2 * Basal_area_2)
    Mod_1 <- sum(density_1 * Basal_area_1)

    # Extract coefficients for the given species
    coef <- coefficients[sp, ]
    growth <- coef$growth
    birth <- coef$birth
    mortality <- coef$mortality
    LC_growth <- coef$LC_growth
    LC_birth <- coef$LC_birth
    LC_mortality <- coef$LC_mortality


    # Calculate new densities
    for (i in 1:length(density_1)) {
        density_1_new <- append(density_1_new,
            max(0, density_1[i] - growth[i] * density_1[i] * (1 - LC_growth[i] * Mod_2) +
            birth[i] * Basal_area_2 * density_2[i] * (1 - LC_birth[i] * (Mod_1 + Mod_2)) -
            density_1[i] * (LC_mortality[i] * Mod_2 + mortality[i]))
        )

      density_2_new <- append(density_2_new,
            max(0, density_2[i] + growth[i] * density_1[i] * (1 - LC_growth[i] * Mod_2) -
            mortality[i] * density_2[i])
        )
    }

    return(list(density_1_new,density_2_new))
}

forest_new_list(c(20,20), c(20,20), c("Resineux","Hetre"))

forest_simul_list <- function(T, density_comp_1, density_comp_2, control, species_selection, CONTROL = TRUE) {
    # forest_simul_list function
    # Input:
    #   T: Number of time steps
    #   density_comp_1: Initial density for compartment 1
    #   density_comp_2: Initial density for compartment 2
    #   control: List of control values
    #   species_selection: Species information
    # Output:
    #   Data frame containing simulation results with columns:
    #   - density_comp_1: Density for compartment 1
    #   - density_comp_2: Density for compartment 2
    #   - Prod: Productivity
    #   - control: Control values
    #   - t: Time step
    #   - Species: Species information

    # Number of species
    N_sp <- length(species_selection)

    # Initialize the forest dataframe
    forest <- data.frame(
        density_comp_1 = density_comp_1,
        density_comp_2 = density_comp_2,
        Prod = 0,
        Control = NA,
        Time = 1,
        Species = species_selection,
        stringsAsFactors = FALSE
    )

    # Main simulation loop
    for (t in 2:T) {
        u <- c(0, 0)
        # Apply control values every 5th time step starting from the 2nd
        if (t %% 5 == 2 & CONTROL) {
            post_cut <- apply_control(t, control, forest, species_selection)
            forest <- forest %>% rbind(post_cut)
        }
        # Generate new values for the forest
        new_forest <- forest_new_list(
                tail(forest$density_comp_1, N_sp), 
                tail(forest$density_comp_2, N_sp),
                species_selection)
        new_forest <- data.frame(new_forest[[1]], new_forest[[2]], 
            0, NA, t, species_selection, stringsAsFactors = FALSE)
        colnames(new_forest) <- colnames(forest)
        forest <- forest %>% rbind(new_forest)
    }

    #forest <- Adjust_output(forest)
    return(forest)
}

apply_control <- function(t, control, forest, species_selection) {
    N_sp <- length(species_selection)
    u <- control[[(t + 3) / 5]]
    actual <- tail(forest$density_comp_2, N_sp)
    # Adjust control values based on the actual values
    u[u > actual] <- actual[u > actual]
    out <- actual - u
    # Create a dataframe for post-cut values
    post_cut <- data.frame(
      as.numeric(tail(forest$density_comp_1, N_sp)),
      u,
      0,
      out,
      t - 1,
      species_selection,
      stringsAsFactors = FALSE
    )
    colnames(post_cut) <- colnames(forest)
    return(post_cut)
}

Adjust_output <- function(forest){
    forest <- forest %>% mutate(is_control = ifelse(is.na(Control), 0, 1))
    short_forest <- forest %>%
        left_join(coefficients %>% select(Species, volume_1, volume_2), by = "Species") %>%
        group_by(is_control, Time, Species) %>%
        mutate(Biomass = density_comp_1 * volume_1 + density_comp_2 * volume_2,
            Prod = volume_2 * Control) %>%
        ungroup() %>%
        group_by(Time, is_control) %>%
        summarize(Biomass = sum(Biomass), Prod = sum(Prod)) %>% ungroup()
    short_forest[which(is.na(short_forest$Prod)),"Prod"] <- 0
    short_forest <- short_forest %>% mutate(Tot_prod = cumsum(Prod))
    species_forest <- forest %>%
        pivot_longer(
            cols = c(density_comp_1, density_comp_2), 
            names_to = "compartiment", 
            values_to = "Density") %>%
        group_by(Time, Species, is_control) %>%
        summarize(Density = sum(Density)) %>%
        group_by(is_control, Time) %>%
        mutate(Density_Tot = sum(Density)) %>%
        mutate(Shannon = Density/Density_Tot * log2(Density/Density_Tot)) %>%
        group_by(Time, is_control) %>%
        summarize(Shannon = -sum(Shannon))
    forest <- forest %>%
        select(-Prod) %>%
        left_join(short_forest, by = c("Time", "is_control")) %>%
        left_join(species_forest, by = c("Time", "is_control"))
    forest <- forest %>% pivot_longer(cols = c(density_comp_1, density_comp_2), names_to = "Compartment", values_to = "Density")
    return(forest)
}


multi_sim <- function(EI_min = c(20,20), EI_max = c(100,100)){
    # Simulate a lot of forest with different mixture and initial states
    all_sp <- list(c("Resineux","Bouleau"), c("Resineux","Hetre"), c("Resineux","Chene"), c("Bouleau","Hetre"),
        c("Bouleau","Chene"), c("Bouleau","Pin"), c("Hetre","Chene"), c("Hetre","Pin"), c("Chene","Pin"))
    all_EI <- list(EI_min, EI_max)
    
    forest <- data.frame()
    for (i in 1:9){
        for(j in 1:2){
        forest <- forest %>% 
            rbind(
                cbind(forest_simul_list(200, all_EI[[j]], all_EI[[j]],u,all_sp[[i]], FALSE), paste(all_sp[[i]][1],"-", all_sp[[i]][2]), j))
        }
    }
    forest <- forest %>% select(- Control, -Prod)
    forest <- forest %>%
        pivot_longer(
            cols = c(density_comp_1, density_comp_2), 
            names_to = "compartiment", 
            values_to = "Density")
    print(head(forest))
    colnames(forest) <- c("Time","Species","Association", "EI", "Compartment", "Density")
    
    ggplot(forest) +
        geom_line(aes(x = Time, y = Density, color = Species, linetype = Compartment)) +
        scale_linetype_manual(values = c("dashed", "solid")) +
        theme(legend.position = "none")+
        labs(x = "Temps", y = "Densité") +
        theme_bw() +
        ylim(0,300) +
        facet_grid(EI ~ Association)
}
multi_sim()

################""""""
# Viability

# forest + 5 ans
forest_5 <- function(density_1, density_2, sp, u){
    Year = forest_new_list(density_1, u, sp)
    for(i in 1:4){
        Year = forest_new_list(Year[[1]], Year[[2]] , sp)
    }
    return(Year)
}

forest_transition_table <- function(species_selection){
    table = expand.grid(x1,x1,x2,x2,u,u)
    colnames(table) = c("x1_1","x1_2","x2_1","x2_2","u1","u2")
    
    table <- table %>%
        filter(u1 <= x2_1 & u2 <= x2_2)

    metric <- table %>% mutate(
        Biomass = coefficients[species_selection[1],]$volume_1 * x1_1 + coefficients[species_selection[2],]$volume_1 * x1_2 +
            coefficients[species_selection[1],]$volume_2 * x2_1 + coefficients[species_selection[2],]$volume_2 * x2_2,
        Shannon = -((x1_1 + x2_1)/(x1_1 + x2_1 + x1_2 + x2_2) * log2((x1_1 + x2_1)/(x1_1 + x2_1 + x1_2 + x2_2)) +
            (x1_2 + x2_2)/(x1_1 + x2_1 + x1_2 + x2_2) * log2((x1_2 + x2_2)/(x1_1 + x2_1 + x1_2 + x2_2))),
        Prod = coefficients[species_selection[1],]$volume_2 * (x2_1 - u1) + coefficients[species_selection[2],]$volume_2 * (x2_2 - u2)) %>%
        select(Biomass, Shannon, Prod)
    
    # change NAn in shannon by 0
    metric[is.na(metric)] <- 0
    
    table <- table %>%
        rowwise() %>%
        mutate(x1 = list(c(x1_1, x1_2)),
            x2 = list(c(x2_1, x2_2)),
            u = list(c(u1, u2))) %>%
        select(x1, x2, u)


    table$x1_new <- apply(table, 1, function(row) {
      result <- forest_5(row[["x1"]], row[["x2"]], species_selection, row[["u"]])
      return(result[1])
    })
    table$x2_new <- apply(table, 1, function(row) {
      result <- forest_5(row[["x1"]], row[["x2"]], species_selection, row[["u"]])
      return(result[2])
    })

    table <- table %>% 
        rowwise() %>%
        mutate(x1_new = list(as.double(unlist(x1_new))),
            x2_new = list(as.double(unlist(x2_new))))
    table <- table %>% cbind(metric)

    states <- Map(c, table$x1, table$x2)
    new_states <- Map(c, table$x1_new, table$x2_new)
    mat_states <- do.call(rbind, states)
    mat_new_states <- do.call(rbind, new_states)
    nn2(mat_states, query = mat_new_states, k = 1)$nn.idx -> idx
    table <- table %>% cbind(table[idx,c("x1","x2")])
    colnames(table) <- c("x1", "x2", "u","x1_new", "x2_new", "Biomass", "Shannon", "Prod", "x1_nn", "x2_nn")
    return(table)
}

fill_parallel <- function(x, species_selection){
  Ncpus <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(cl)
  
  res <- foreach::foreach(i=1:nrow(x), .combine='rbind', .packages = c("tidyverse"), .export = c("forest_5", "forest_new_list", "coefficients")) %dopar%{
          return(forest_5(x[[i, "x1"]][[1]], x[[i, "x2"]][[1]], species_selection, x[[i, "u"]][[1]]))
  }
  
  parallel::stopCluster(cl)
  return(res)
}

fill_parallel(table, species_selection)

fill_transition_matrix <- function(table){
    M = array(0, dim = c(
        length(x1), length(x1), length(x1), length(x1), # départ
        length(x1), length(x1), length(x1), length(x1), # arrivée
        length(u), length(u))) # control

    for(i in 1 : nrow(table)){
        M[table[[i,"x1"]][1]/f + 1, table[[i,"x1"]][2]/f + 1, table[[i,"x2"]][1]/f + 1, table[[i,"x2"]][2]/f +1,
        table[[i,"x1_nn"]][1]/f + 1, table[[i,"x1_nn"]][2]/f + 1, table[[i,"x2_nn"]][1]/f + 1, table[[i,"x2_nn"]][2]/f +1,
            table[[i,"u"]][1]/fc +1, table[[i,"u"]][1]/fc +1] = 1
    }
    return(M)
}

make_viability_kernel <- function(M, table){
    all_control = expand.grid(1:length(u),1:length(u))
    V = array(0, dim = c(length(x1), length(x1), length(x1), length(x1), length(u), length(u)))

    for (i in 1:nrow(table)){
            V[table[[i,"x1"]][1]/f + 1, table[[i,"x1"]][2]/f + 1, table[[i,"x2"]][1]/f + 1, table[[i,"x2"]][2]/f +1,
            table[[i,"u"]][1]/f + 1, table[[i,"u"]][2]/f + 1] = 
                ifelse(table$Shannon[i] > 0.6 & table$Biomass[i] > 10 & table$Prod[i] > 0, 1, 0)
    }

    for(t in 1:50){
        V_new = array(0, dim = dim(V))
        for (i in 1:length(x1)){
            for (j in 1:length(x1)){
                for (k in 1:length(x1)){
                    for (l in 1:length(x1)){
            if(sum(V[i,j,k,l,,] == 0)){next}
            for (v in 1:length(u)){
                for(w in 1:length(u)){
                #v = all_control[u,1]
                #w = all_control[u,2]
                    if(sum(M[i,j,k,l,,,,,v,w] * V[,,,,v,w]) == 1){V_new[i,j,k,l,v,w] = 1; break}
                }}}}}
        }
        if(all(V == V_new)){print(t);break}
        V = V_new
    }
    return(V)
}


# Ajouter la viabilité
table_viability <- function(table, sh = 0.6, biom = 10, prod = 30){
    # choix de contraintes
    table <- table %>% mutate(Viable = ifelse(Shannon > sh & Biomass > biom & Prod > prod, 1, 0), Constraints = Viable)
    diff = TRUE
    i = 1
    while(diff & i < 20){
        i = i+1
        Viability <- table %>% select(x1,x2,Viable) %>% group_by(x1,x2) %>% summarize(Viable = sum(Viable)) %>%
            mutate(Viable = ifelse(Viable == 0, 0, 1))
        colnames(Viability) <- c("x1_nn", "x2_nn", "Viable_nn")
        table <- table %>% left_join(Viability, by = c("x1_nn","x2_nn"))
        table <- table %>% mutate(Viable_new = Viable * Viable_nn)
        table <- table %>% select(-Viable_nn)
        diff = !all(table$Viable == table$Viable_new)
        table$Viable = table$Viable_new
    }
    print(paste("calculated whith : ", i, "iterations"))
    return(table)
}

# Visualisons une dimension de V
#V1 <- apply(V, c(1,2,3,4), sum)
#V1[V1 > 0] <- 1
#V1 <- V1[1,1,,]
#image(V1)
