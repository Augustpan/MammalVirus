library(tidyverse)
library(vegan)
library(ggsci)
library(scatterpie)
library(ggConvexHull)

getPathToVirusList = function() { "data/mammal_virus_list.csv" }
getPathToHostList = function() { "data/host_occurence.csv" }
getPathToHostInfo = function() { "data/host_metadata.csv" }

getPathToAnnotatedVirusList = function() { "data/mammal_virus_list_annotated.csv" }

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