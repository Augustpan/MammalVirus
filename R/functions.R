library(tidyverse)
library(vegan)
library(ggsci)

getPathToVirusList = function() { "data/mammal_virus_list.csv" }
getPathToHostList = function() { "data/host_occurence.csv" }
getPathToHostInfo = function() { "data/host_metadata.csv" }

# merge host_list and host_info tables and transform into site-by-species table
buildHostCommunityTable = function(host_list, host_info) {
  host_merged = inner_join(host_info, host_list, by = "Mammal_species")
  host_comm = host_merged %>%
    select(Jingmen, Longquan, Wenzhou, Wufeng) %>%
    t()
  colnames(host_comm) = host_merged$Mammal_species
  host_comm
}

calcHostAlphaDiversity = function(host_comm, host_info, host, method, nsample=1) {
  if (host == "all") { 
    host_filter = host_info$Mammal_species 
  } else { 
    filter(host_info, Mammal==host)$Mammal_species
  }

  sample_specnum = matrix(nrow=nsample, ncol=4)
  for (i in 1:nsample) {
    tmp = host_comm %>%
      rrarefy(620) %>%
      as_tibble() %>%
      select(all_of(host_filter)) %>%
      method()
    sample_specnum[i,] = tmp
  }
  colnames(sample_specnum) = rownames(host_comm)
  
  tb_specnum = sample_specnum %>%
    t() %>%
    as.data.frame() %>%
    mutate(site_name = rownames(.)) %>%
    pivot_longer(cols = 1:nsample, names_to="nsample")
}

createPlotHostAlphaDiversity = function(tb_specnum) {
  ggplot(tb_specnum) + 
    geom_bar(aes(x=site_name, y=value, fill=site_name), stat = "identity") +
    scale_fill_uchicago() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
    ylab("Host Richness") + 
    xlab("")
}

# CA plot of viral organophilism
createPlotCaOrganophilism = function(virus_list) {
  virus_list %>%
    select(c("lung_Abundance",
             "liver_Abundance",
             "spleen_Abundance",
             "kindey_Abundance",
             "feces_Abundance")) %>%
    vegan::cca() %>%
    plot()
}