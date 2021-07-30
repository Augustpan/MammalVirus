library(tidyverse)


#### MAKE CLEAN VIRUS LIST

seg = read_tsv("data/segment.tsv", show_col_types = FALSE)
virus_list = read_csv("data/mammal_virus_list.csv", show_col_types = FALSE)
sp = str_split(virus_list$Host_Name, " ", simplify = TRUE)
virus_list$Site_Name = sp[,1]
virus_list$Host_Species = paste0(sp[,2], "_", sp[,3])
virus_list$Virus_Species = virus_list$Viral_Family_Genus_Species

multi_host_virus = virus_list %>%
    select(Virus_Species, Host_Species) %>%
    distinct() %>%
    count(Virus_Species) %>%
    filter(n > 1) %>%
    `$`(Virus_Species)

virus_list %>%
    select(
        Site_Name, 
        Viral_Family, Viral_Genus, Virus_Species, 
        Host_Order=Mammal, Host_Species,
        Genome_Type = Genome_type,
        Novelty, Zoonosis,
        Lung_Abundance = lung_Abundance,
        Kidney_Abundance = kindey_Abundance,
        Liver_Abundance = liver_Abundance,
        Feces_Abundance = feces_Abundance,
        Spleen_Abundance = spleen_Abundance,
        Total_Abundance) %>%
    left_join(seg, by="Viral_Family") %>%
    mutate(Multi_Host = Virus_Species %in% multi_host_virus) %>%
    #mutate_if(is.character, as.factor) %>%
    saveRDS(file="data/virus_list.rds")

#### MAKE VIRAL FAMILY DISTANCE MATRIX

virus_list %>%
    select(Viral_Family) %>%
    distinct() %>%
    `$`(Viral_Family)

#### MAKE HOST ORDER DISTANCE MATRIX

host_tree = read.tree("data/host_order.nwk")
host_phylo_corr_matrix = vcv.phylo(host_tree, corr=TRUE)
saveRDS(host_phylo_corr_matrix, file="data/hostOrdPrecisionMat.rds")