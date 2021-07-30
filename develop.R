seq_stat = tb_seq_stat
virus_list = tb_virus_list_clean %>%
    pivot_longer(
        cols=Lung_Abundance:Spleen_Abundance, 
        names_to="Sample_type", 
        values_to = "Abundance") %>%
    mutate(Sample_type = str_split(Sample_type, "_", simplify = TRUE)[,1] %>% str_to_lower())

sp = seq_stat %>%
    `$`(`Host_species (No. of individual)`) %>%
    str_split(" ", simplify = TRUE)

seq_stat$Host_Species = paste0(sp[,1], "_", sp[,2])

full_list = inner_join(virus_list, seq_stat, by=c("Host_Species", "Site_Name"="Sampling_city", "Sample_type"))

full_list$reads = round(full_list$Abundance * full_list$`norRNA_data (total_read_counts)`)

pool = full_list %>%
    select(Site_Name, Host_Species, Sample_Type=Sample_type, Virus_Species, Reads=reads) %>%
    pivot_wider(names_from=Virus_Species, values_from=Reads, values_fill=0)

host = pool %>%select(-Sample_Type) %>%group_by(Site_Name, Host_Species) %>% summarise_all(sum) %>% ungroup()
site = host %>%select(-Host_Species) %>% group_by(Site_Name) %>% summarise_all(sum) %>% as.data.frame()

rownames(site) = site$Site_Name


df = virus_list %>% 
    group_by(Site_Name, Host_Species, Virus_Species, Multi_Host) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(Site_Name, Host_Species) %>%
    summarise(Multi_Host = max(Multi_Host) == 1, Richness=length(unique(Virus_Species))) %>%
    inner_join(
        virus_list %>% 
            group_by(Site_Name, Host_Species, Virus_Species, Multi_Host) %>%
            summarise(Abundance = sum(Abundance)) %>%
            ungroup() %>%
            group_by(Site_Name) %>%
            summarise(Regional_Richness=length(unique(Virus_Species))),
        by = "Site_Name"
    )
