
# Function for the viability
#--------------------------------------
########################################

transition_table <- function(N_layers = 3, species_selection, N_sp = 2, N_cut = 2, max_tree = c(800, 400, 200), pas = 2){

    Grid_state <- state_combinatory(N_layers, species_selection, N_sp, N_cut, max_tree, pas)
    Metric = metric_parallel(Grid_state, N_layers, species_selection)
    Final_state <- fill_parallel(Grid_state, N_layers, species_selection)
    
    Init_state_idx <- 1:nrow(Grid_state)
    Final_state_idx <- nn2(Grid_state, Final_state, k = 1)$nn.idx

    Control_state <- cbind(Grid_state, Final_state_idx, 
        control_vol_wood = Metric$vol_wood,
        control_shannon = Metric$shannon,
        control_shannon_vert = Metric$shannon_vert,
        control_density_tot = Metric$density_tot,
        control_basal_area = Metric$basal_area)
    Init_state <- cbind(Grid_state, Init_state_idx, Metric)

    table <- merge(Init_state, Control_state, by = c(1:(N_sp * (N_layers - N_cut))))

    column_filter = paste0("Var", (N_sp * (N_layers - N_cut) + 1): (N_sp * (N_layers - N_cut) + N_sp * N_cut))

    table <- table[rowSums(table[paste0(column_filter, ".x")] >= table[paste0(column_filter, ".y")]) == length(column_filter), ]
    table <- table %>% mutate(extraction = vol_wood - control_vol_wood)

    return(table)
}

state_combinatory <- function(N_layers = 3, species_selection, N_sp = 2, N_cut = 2,  max_tree = c(800, 400, 200), pas = 2){
    init1 <- seq(0, max_tree[1], max_tree[1]/pas)
    init2 <- seq(0, max_tree[2], max_tree[2]/pas)
    init3 <- seq(0, max_tree[3], max_tree[3]/pas)

    table <- expand.grid(c(rep(list(init1), N_sp), rep(list(init2), N_sp), 
    rep(list(init3), N_sp)))
    return(table)
}

columns <- function(State, N_layers, N_cut, N_sp){
    if(State == "Control_State"){
        return(c(1:((N_layers-N_cut)*N_sp), (N_layers * N_sp + 1): (N_layers * N_sp + N_sp * N_cut)))
    }
    if(State == "Initial_State"){
        return(c(1:(N_layers * N_sp)))
    }
    if(State == "Final_State"){
        return(c((N_layers * N_sp + N_sp * N_cut + 1): (N_layers * N_sp + N_sp * N_cut + N_sp * N_layers)))
    }
}

fill_parallel <- function(x, N_layers, species_selection){
    Ncpus <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(Ncpus)
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=1:nrow(x), .combine='rbind', .packages = c("tidyverse"),
      .export = c("forest_new_5", "forest_new", "forest_new_table", "coefficients", "basal_area")) %dopar%{
        return(forest_new_table(coefficients, x[i,], N_layers, species_selection))
    }
  parallel::stopCluster(cl)
  return(res)
}

forest_new_table <- function(coefficients,row_table, N_layers, species_selection){
    N_sp = length(species_selection)
    density = matrix(row_table %>% as.numeric(), nrow = N_layers, ncol = N_sp, byrow = TRUE) %>% data.frame()
    colnames(density) = species_selection
    density = forest_new_5(coefficients,density) %>% data.frame() %>% unlist() %>% as.numeric()
    return(density)
}

metric_parallel <- function(x, N_layers, species_selection){
    Ncpus <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(Ncpus)
    doParallel::registerDoParallel(cl)
    res <- foreach::foreach(i=1:nrow(x), .combine='rbind', .packages = c("tidyverse"),
      .export = c("forest_metric_table", "forest_metric", "basal_area", "Volume")) %dopar%{
        return(forest_metric_table(x[i,], N_layers, species_selection))
    }
  parallel::stopCluster(cl)
  colnames(res) <- c("density_tot", "basal_area", "vol_wood", "shannon_vert", "shannon")
  res[is.na(res)] <- 0
  return(data.frame(res))
}

forest_metric_table <- function(row_table, N_layers, species_selection){
    N_sp = length(species_selection)
    density = matrix(row_table %>% as.numeric(), nrow = N_layers, ncol = N_sp, byrow = TRUE) %>% data.frame()
    colnames(density) = species_selection
    metric = forest_metric(density)
    return(metric)
}

table_viability <- function(table,
    vol_wood =0, ba =0, shannon =0, ext =0, shannon_vert =0,
    control_vol_wood =0, control_shannon =0, control_shannon_vert =0,
    control_density_tot =0, control_basal_area = 0){
    # choix de contraintes
    table$Viable <- as.numeric(
        table$vol_wood >= vol_wood &
        table$basal_area >= ba &
        table$shannon >= shannon &
        table$extraction >= ext &
        table$shannon_vert >= shannon_vert,
        table$control_vol_wood >= control_vol_wood &
        table$control_shannon >= control_shannon &
        table$control_shannon_vert >= control_shannon_vert &
        table$control_density_tot >= control_density_tot &
        table$control_basal_area >= control_basal_area)
    
    V <- table %>% select(Init_state_idx, Final_state_idx, Viable)
    V$Constraints  = V$Viable
    for(i in 1:20){
        Viability <- V %>% select(Init_state_idx,Viable) %>% group_by(Init_state_idx) %>% summarize(Viable = sum(Viable)) %>%
            mutate(Viable = ifelse(Viable == 0, 0, 1))
        colnames(Viability) <- c("Final_state_idx", "Viable_nn")
        V <- V %>% left_join(Viability, by = c("Final_state_idx"))
        V <- V %>% mutate(Viable_new = Viable * Viable_nn)
        if(all(V$Viable == V$Viable_new)){break} 
        V <- V %>% select(-Viable_nn)
        V$Viable = V$Viable_new
    }
    print(paste("calculated whith : ", i, "iterations"))
    return(V %>% select(Constraints, Viable, Viable_nn))
}


LCg = 0.005903714
LCb = 0.006468877
LCm = 0.004589390
coefficients <- data.frame(
    growth1 = c(0.5, 0.5, 0.29, 0.5, 0.37, 0.28),
    growth2 = c(0.96, 0.9, 0.53, 0, 0.71, 0.5),
    birth = c(60, 60, 60, 60, 6, 60),
    mortality = c(0.013, 0.015, 0.023, 0.031, 0.012, 0.008),
    LC_growth = LCg * c(1, 5, 9, 9, 1, 7),
    LC_birth = LCb* c(3, 5, 7, 7, 3, 7),
    LC_mortality = LCm * c(1, 5, 9, 9, 1, 7)
)
rownames(coefficients) <- c("Sapin", "Epicea", "Pin", "Bouleau", "Hetre", "Chene")
colnames(coefficients) <- c("growth1","growth2", "birth", "mortality", "LC_growth", "LC_birth", "LC_mortality")

basal_area = c(0.01, 0.16, 0.9)

Volume = Volume = readRDS("Data/volume_table.rds")

latin = c("Abies alba", "Picea abies", "Pinus sylvestris", "Betula pendula", "Fagus sylvatica", "Quercus pubescens")
common = c("Sapin", "Epicea", "Pin", "Bouleau", "Hetre", "Chene")
contraction = c("AAlb", "Pabi", "PSyl", "BPen", "FSyl", "QPub")
latin_name <- data.frame(latin, common, contraction)
rownames(latin_name) <- common