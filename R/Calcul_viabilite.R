source("Fonctions.R")
library(plotly)


table <- transition_table(2, c("Resineux","Bouleau"), 1, n_reprod = 1, max_tree = 200, freq = 5)
saveRDS(table, "Data/table_27_10.RDS")

Shannon = seq(0,1,0.05)
Biomass = seq(0,80,2)
Extraction = seq(0,40,2)
sensibilite <- data.frame()
# tester la sensibilité aux différentes contraintes
for(s in Shannon){
    for(b in Biomass){
        for(p in Extraction){
            table_v <- table_viability(table, b, s, p, 0)
            sensibilite <- rbind(sensibilite,
                c(s,b,p,table_v %>% filter(Constraints == 1) %>% nrow(),
                table_v %>% filter(Viable == 1) %>% nrow()))
        }
    }
}
colnames(sensibilite) <- c("Shannon", "Biomass", "Extraction", "Constraints", "Viable")
saveRDS(sensibilite, "sensibilite_long.RDS")