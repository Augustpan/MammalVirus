library(tidyverse)
library(vegan)
library(ggsci)
library(scatterpie)
library(ggConvexHull)
library(UpSetR)
library(ape)
library(Matrix)

getPathToVirusList = function() { "data/mammal_virus_list.csv" }
getPathToHostList = function() { "data/host_occurence.csv" }
getPathToHostInfo = function() { "data/host_metadata.csv" }
getPathToHostSpeciesTree = function() { "data/host_phylo_tree.nwk" }
getPathToHostOrderTree = function() { "data/host_order.nwk" }
getPathToVirusSegment = function() { "data/segment.tsv" }
getPathToVirusFamPreMat = function() { "data/virusFamPrecisionMat.rds" }
getPathToHostOrdPreMat = function() { "data/hostOrdPrecisionMat.rds" }
getPathToSeqStat = function() { "data/lib_seq_stat.csv" }

cleanVirusList = function(virus_list, virus_segment) {
  sp = str_split(virus_list$Host_Name, " ", simplify = TRUE)
  
  # standardize some column names
  virus_list$Site_Name = sp[,1]
  virus_list$Host_Species = paste0(sp[,2], "_", sp[,3])
  virus_list$Virus_Species = virus_list$Viral_Family_Genus_Species
  
  # find out the viruses that were found in multiple host species
  multi_host_virus = virus_list %>%
    select(Virus_Species, Host_Species) %>%
    distinct() %>%
    count(Virus_Species) %>%
    filter(n > 1) %>%
    `$`(Virus_Species)
  
  # select and rename columns, joining segmentation info
  virus_list = virus_list %>%
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
    left_join(virus_segment, by="Viral_Family") %>%
    mutate(Multi_Host = Virus_Species %in% multi_host_virus)
  
  # rename Order of host species
  f1 = virus_list$Host_Order == "shrews"
  f2 = virus_list$Host_Order == "bats"
  f3 = virus_list$Host_Order == "rodents"
  virus_list$Host_Order[f1] = "Insectivora"
  virus_list$Host_Order[f2] = "Chiroptera"
  virus_list$Host_Order[f3] = "Rodentia"

  # drop viral Family that there is not enough replications
  to_drop = c("Adenoviridae", "Anelloviridae", "Arteriviridae", "Circoviridae", "Matonaviridae")
  virus_list = filter(virus_list, !(Viral_Family %in% to_drop))
  
  # viral richness within host population
  pop_viral_richness = virus_list %>% count(Site_Name, Host_Order, Host_Species, name="Population_Viral_Richness")
  
  # viral richness within a sample site (or a host community)
  com_viral_richness = virus_list %>% 
    count(Site_Name, Host_Species, Virus_Species) %>%
    select(-Host_Species) %>%
    distinct() %>%
    count(Site_Name, name="Region_Viral_Richness")
  
  # host richness within a sample site
  com_host_richness = tibble(
    Site_Name = c("Longquan","Wufeng","Wenzhou","Jingmen"),
    Region_Host_Richness = c(17,18,17,13)
  )
  
  # merge richness info in to virus_list
  virus_list = virus_list %>%
    inner_join(pop_viral_richness, by=c("Site_Name", "Host_Order", "Host_Species")) %>%
    inner_join(com_viral_richness, by="Site_Name") %>%
    inner_join(com_host_richness, by="Site_Name")
  
  # virus tissue tropism
  ttca = virus_list %>% select(10:14) %>% vegan::cca()
  virus_list$TT_CA1 = ttca$CA$u[,1]
  virus_list$TT_CA2 = ttca$CA$u[,2]
  
  # viral community composition within host population
  tmp = virus_list %>% 
    pivot_wider(names_from = Viral_Genus, values_from=Total_Abundance, values_fill=0) %>% 
    select(Site_Name, Host_Species,21:ncol(.)) %>% 
    group_by(Site_Name, Host_Species) %>% 
    summarise_all(sum) %>%
    ungroup()
  vgca = tmp %>%
    select(3:ncol(.)) %>%
    vegan::cca()
  virus_list = virus_list %>% inner_join(
    tibble(
      VG_CA1 = vgca$CA$u[,1],
      VG_CA2 = vgca$CA$u[,2],
      Site_Name = tmp$Site_Name,
      Host_Species = tmp$Host_Species),
    by = c("Site_Name", "Host_Species"))
  
  # transform strings into factors
  virus_list = mutate_if(virus_list, is.character, as.factor)
  
  # return
  virus_list
}

makeFullModelList = function(virus_list, virusFamPrecisionMat, hostOrdPrecisionMat) {
  virusFamPrecisionMat = virusFamPrecisionMat %>%
    `[`(unique(virus_list$Viral_Family), unique(virus_list$Viral_Family)) %>%
    solve()
  virusFamRank = virusFamPrecisionMat %>%
    rankMatrix() %>%
    as.vector()
  hostOrdPrecisionMat = hostOrdPrecisionMat %>%
    `[`(c("Insectivora", "Chiroptera", "Rodentia"), c("Insectivora", "Chiroptera", "Rodentia")) %>%
    solve()
  hostOrdRank = hostOrdPrecisionMat %>%
    rankMatrix() %>%
    as.vector()
  
  # the biological characteristics of virus
  effect_viral_family = "s(Viral_Family, bs = 're', xt = list(S = list(virusFamPrecisionMat), rank = virusFamRank))"
  effect_tissue_tropism_axis1 = "s(TT_CA1, bs = 'tp', k = 10)"
  effect_tissue_tropism_axis2 = "s(TT_CA2, bs = 'tp', k = 10)"
  effect_genome_type = "Genome_Type"
  effect_zoonosis = "Zoonosis"
  effect_novelty = "Novelty"
  
  # the ecological characteristics of virus
  effect_com_viral_richness = "s(Region_Viral_Richness, bs = 'tp', k = 3)"
  effect_pop_viral_richness = "s(Population_Viral_Richness, bs = 'tp', k = 9)"
  effect_viral_com_genus_axis1 = "s(VG_CA1, bs = 'tp', k = 10)"
  effect_viral_com_genus_axis2 = "s(VG_CA2, bs = 'tp', k = 10)"
  
  # the characteristics of host
  effect_host_order = "s(Host_Order, bs = 're')"
  effect_host_order_specific_viral_family = "s(Host_Order, Viral_Family, bs = 're')"
  effect_host_order_phylogeny = "s(Host_Order, bs = 're', xt = list(S = list(hostOrdPrecisionMat), rank = hostOrdRank))"
  effect_com_host_richness = "s(Region_Host_Richness, bs = 'tp', k = 3)"
  
  # model terms
  Terms = list(
    HostEffects = c(NA,
                    effect_host_order,
                    effect_host_order_specific_viral_family,
                    effect_host_order_phylogeny),
    TissueTropismEffects1 = c(NA, effect_tissue_tropism_axis1),
    TissueTropismEffects2 = c(NA, effect_tissue_tropism_axis2),
    CommunityCompositionEffects1 = c(NA, effect_viral_com_genus_axis1),
    CommunityCompositionEffects2 = c(NA, effect_viral_com_genus_axis2),
    ViralFamilyEffect = c(NA, effect_viral_family), 
    GenomeTypeEffect = c(NA, effect_genome_type),
    ZoonosisEffect = c(NA, effect_zoonosis),
    NoveltyEffect = c(NA, effect_novelty),
    
    RegionalViralRichnessEffect = c(NA, effect_com_viral_richness),
    PopulationViralRichnessEffect = c(NA, effect_pop_viral_richness),
    HostRichnessEffect = c(NA, effect_com_host_richness)
  )
  
  CompetingFullModels = expand.grid(Terms) %>% 
    filter(!(HostEffects == "s(Host_Order, Viral_Family, bs = 're')" & is.na(ViralFamilyEffect))) %>%
    apply(1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
    data.frame(Formula = ., stringsAsFactors = F) %>% 
    mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
           Formula = paste('Multi_Host ~', Formula))
  
  # return
  CompetingFullModels
}

fitModel = function(formula, data) {
  fit = gam(
    as.formula(formula),
    data = data,
    family = "binomial",
    method = "REML",
    select = FALSE
  )
  fit
}

seqDepthValidation = function() {
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
}

calcHostAlphaDiversity = function(host_list, host_info, method=specnumber, nsample=1) {
  host_merged = inner_join(host_info, host_list, by = "Mammal_species")
  host_comm = host_merged %>%
    select(Jingmen, Longquan, Wenzhou, Wufeng) %>%
    t()
  colnames(host_comm) = host_merged$Mammal_species
  
  all_hosts = host_info$Mammal_species
  rodents = filter(host_info, Mammal=="rodents")$Mammal_species
  bats = filter(host_info, Mammal=="bats")$Mammal_species
  shrews = filter(host_info, Mammal=="shrews")$Mammal_species
  lst_host = list("all"=all_hosts, "rodents"=rodents, "bats"=bats, "shrews"=shrews)
  lst_result = list()
  for (j in 1:length(lst_host)) {
    host_filter = lst_host[[j]]
    host_name = names(lst_host)[j]
    sample_specnum = matrix(nrow=nsample, ncol=4)
    for (i in 1:nsample) {
      tmp = host_comm %>%
        as_tibble() %>%
        select(all_of(host_filter)) %>%
        method()
      sample_specnum[i,] = tmp
    }
    colnames(sample_specnum) = rownames(host_comm)
    rownames(sample_specnum) = 1:nsample
    tb_specnum = sample_specnum %>%
      t() %>%
      as.data.frame() %>%
      mutate(site_name = rownames(.)) %>%
      pivot_longer(cols = 1:nsample, names_to="sample", values_to="richness") %>%
      mutate(species = host_name)
    lst_result[[j]] = tb_specnum
  }
  bind_rows(lst_result)
}

sortVirusTableCA = function(virus_list) {

  tmp = virus_list %>%
    select(Host_Name, Mammal, Viral_Family_Genus, Total_Abundance) %>%
    mutate(Incidence = Total_Abundance > 0) %>%
    group_by(Host_Name, Viral_Family_Genus) %>%
    summarise(Incidence = sum(Incidence) >= 1) 
  vtable = tmp %>%
    left_join(select(virus_list, Host_Name, Mammal) %>% distinct(), by="Host_Name") %>%
    pivot_wider(names_from=Viral_Family_Genus, values_from=Incidence, values_fill=0) %>%
    ungroup()
  
  ca_vc = cca(select(vtable, -Host_Name, -Mammal))
  site_score = ca_vc$CA$u[,1]
  species_score = ca_vc$CA$v[,1]
  
  dft = select(vtable, -Host_Name, -Mammal) %>%
    as.data.frame()
  dft[order(site_score),order(species_score)] %>%
    cbind(vtable$Host_Name[order(site_score)], vtable$Mammal[order(site_score)])
}

createPlotSharedHostSpecies = function(host_list, host_info) {
  host_merged = inner_join(host_info, host_list, by = "Mammal_species")
  host_merged = mutate(host_merged, SH=(Jingmen > 0) + (Longquan > 0) + (Wenzhou > 0) + (Wufeng > 0))
  SH = host_merged %>%
    filter(SH > 1)
  SH$x = 1:nrow(SH)
  SH$tot = SH$Jingmen + SH$Wufeng + SH$Wenzhou + SH$Longquan
  ggplot() + 
    geom_scatterpie(
      aes(
        x=x*15, 
        y=0, 
        group=Mammal_species, 
        r = log(tot)
        ), 
      data=SH,
      cols=c("Jingmen", "Wufeng", "Wenzhou", "Longquan")) + 
     geom_text(aes(x=x*15+15, y = 45, label=Mammal_species), data=SH, angle=45) + 
    coord_equal()
}

createPlotCAVirusOrganophilism = function(virus_list) {
  ca_op = virus_list %>%
    select(c("lung_Abundance",
             "liver_Abundance",
             "spleen_Abundance",
             "kindey_Abundance",
             "feces_Abundance")) %>%
    cca()
  
  # species_score
  u = ca_op$CA$u
  
  # site_score
  v = ca_op$CA$v
  lb = c("Lung", "Liver", "Spleen", "Kidney", "Feces")
  lst = c("Arena_Mammarenavirus@Wenzhou mammarenavirus",
    "Paramyxo_Henipavirus@Novel henipavirus 1",
    "Flavi_Hepacivirus@Novel hepacivirus 3",
    "Hanta_Orthohantavirus@Seoul orthohantavirus",
    "Hanta_Orthohantavirus@Hantaan orthohantavirus",
    "Paramyxo_Jeilongvirus@Novel jeilongvirus 3",
    "Paramyxo_Jeilongvirus@Beilong jeilongvirus",
    "Flavi_Hepacivirus@Novel hepacivirus 2",
    "Flavi_Pegivirus@Novel pegivirus 2",
    "Flavi_Pegivirus@Novel pegivirus 3",
    "Flavi_Pegivirus@Novel pegivirus 1")
  dft = filter(virus_list, Viral_Family_Genus_Species %in% lst)
  ggplot() +
    geom_point(aes(x=u[,1], y=u[,2], color=virus_list$Viral_Family), size=3) +
    geom_point(aes(x=u[dft$`No.`,1], y=u[dft$`No.`,2]), shape=2, size=5) +
    #stat_ellipse(aes(x=u[,1], y=u[,2] ,color=virus_list$Viral_Family)) + 
    #geom_convexhull(aes(x=u[,1], y=u[,2], color=virus_list$Viral_Family, fill=NULL), alpha=0) + 
    geom_segment(aes(x=0, y=0, xend=v[,1], yend=v[,2]), arrow=arrow()) +
    geom_text(aes(x=v[,1], y=v[,2], label=lb)) +
    xlab("CA Axis1 (44.73%)") + 
    ylab("CA Axis2 (25.36%)") +
    theme_bw() + 
    theme(axis.title = element_text(size=14, face="bold"),
          axis.text = element_text(size=12),
          legend.text = element_text(size=12))
  
  
  lst_op = list(
    lung = filter(virus_list, lung_Abundance > 0)$Viral_Family_Genus_Species,
    liver = filter(virus_list, liver_Abundance > 0)$Viral_Family_Genus_Species,
    spleen = filter(virus_list, spleen_Abundance > 0)$Viral_Family_Genus_Species,
    kidney = filter(virus_list, kindey_Abundance > 0)$Viral_Family_Genus_Species,
    feces = filter(virus_list, feces_Abundance > 0)$Viral_Family_Genus_Species
  )
  upset(fromList(lst_op))

}

performPERMANOVAVirusOp = function(virus_list) {
  col_names = c("lung_Abundance",
                "liver_Abundance",
                "spleen_Abundance",
                "kindey_Abundance",
                "feces_Abundance")
  vtable = virus_list %>%
    select(all_of(col_names), Viral_Family_Genus_Species, Host_Name) %>%
    pivot_longer(cols=col_names, names_to="organ") %>%
    pivot_wider(names_from=Viral_Family_Genus_Species, values_from=value, values_fill=0)
  Y = vtable[3:ncol(vtable)]
  perm <- how(nperm = 999, blocks = vtable$Host_Name)
  adonis2(Y~organ, data=vtable, method="jaccard", permutations = perm)
  
}

## rewrite needed
virusOP = function() {
  tmp_vir = tb_vir %>%
    select(Viral_Family_Genus_Species, Host_Name, ends_with("_Abundance"), -Max_Abundance, -Total_Abundance)
  
  tmp_lib = tb_lib %>%
    select(
      site_name = Sampling_city,
      host_species = `Host_species (No. of individual)`, 
      reads = `norRNA_data (total_read_counts)`,
      sample_type = Sample_type)
  
  tmp = str_split(tmp_lib$host_species, " ", simplify = TRUE)
  sp_name = paste0(tmp_lib$site_name, " ", tmp[,1], " " ,tmp[,2])
  tmp_lib$Host_Name = sp_name
  tmp_vir = pivot_longer(tmp_vir, cols=3:7, names_to="sample_type", values_to = "Abundance")
  tmp_vir$sample_type = str_split(tmp_vir$sample_type, "_", simplify = TRUE)[,1]
  
  tmp = inner_join(tmp_lib, tmp_vir, by=c("Host_Name","sample_type")) %>%
    mutate(value = value*reads) %>%
    group_by(site_name, Viral_Family_Genus_Species) %>%
    summarise(value=sum(value), reads=sum(reads)) %>%
    pivot_wider(names_from = Viral_Family_Genus_Species, values_from = value, values_fill=0)
  
  tmp_table = tmp %>%
    ungroup() %>%
    select(-all_of(c("site_name","reads"))) %>% 
    round()
}

summariseVirusAnnotation = function(virus_list) {
  list(
    novelty = virus_list %>%
      group_by(Novelty) %>%
      summarise(nov = n()),
    genome = virus_list %>%
      group_by(Genome_type) %>%
      summarise(gen = n()),
    host = virus_list %>%
      group_by(Host_type) %>%
      summarise(hos = n()),
    zoonosis = virus_list %>%
      group_by(Zoonosis) %>%
      summarise(zoo = n())
  )
}

calcHostPhyloDist = function(host_tree, virus_list) {
  host_phylo_dist_mat = cophenetic.phylo(host_tree)
  
  host_phylo_dist = host_phylo_dist_mat %>%
    as.data.frame() %>%
    mutate(A = rownames(.)) %>%
    pivot_longer(cols=1:(ncol(.)-1), names_to="B")
  
  
  sp = str_split(virus_list$Host_Name, " ", simplify = TRUE)
  virus_list$Host_Species = paste0(sp[,2], "_", sp[,3])
  x = filter(virus_list, Viral_Family_Genus_Species %in% lst) %>%
    select(Viral_Family_Genus_Species, Host_Species)
  
  ret_mat = matrix(nrow=length(lst), ncol=6)
  for (i in 1:length(lst)) {
    tmp = filter(x, Viral_Family_Genus_Species == lst[i])$Host_Species
    a = filter(host_phylo_dist, A %in% tmp & B %in% tmp, value > 0)$value %>% summary()
    ret_mat[i,] = a
  }
  colnames(ret_mat) = names(a)
  ret_tb = as_tibble(ret_mat) %>%
    mutate(virus_species = lst)
}

plotVirusGenusComposition = function() {
  ca = sub_vl %>% 
    group_by(Site_Name, Host_Species, Multi_Host, Viral_Genus) %>% 
    summarise(Total_Abundance = sum(Total_Abundance))  %>%
    pivot_wider(names_from=Viral_Genus, values_from=Total_Abundance, values_fill=0) %>%
    group_by(Site_Name, Host_Species) %>% 
    summarise(across(4:ncol(.)-2, ~ sum(.x)), Multi_Host=max(Multi_Host)) %>%
    ungroup() %>%
    select(3:50) %>%
    vegan::cca()
  
  k = sub_vl %>% 
    group_by(Site_Name, Host_Species, Multi_Host, Viral_Genus) %>% 
    summarise(Total_Abundance = sum(Total_Abundance))  %>%
    pivot_wider(names_from=Viral_Genus, values_from=Total_Abundance, values_fill=0) %>%
    group_by(Site_Name, Host_Species) %>% 
    summarise(across(4:ncol(.)-2, ~ sum(.x)), Multi_Host=max(Multi_Host)) %>%
    ungroup() 
  
  mm =c("rodents", "bats", "bats", "rodents", "bats", 
        rep("rodents", 4), "bats", "bats", "rodents",
        "bats", "bats", "rodents", "rodents", "shrews", 
        "shrews", rep("rodents", 3), "bats", "bats")
  
  ggplot() +
    geom_point(aes(x=x, y=y, color=as.factor(k$Multi_Host), shape=mm), size=3) +
    ylim(c(-1,1.5)) +
    xlab("CA Axis1 (11.62%)") + 
    ylab("CA Axis (10.35%)") + 
    theme_bw()
}