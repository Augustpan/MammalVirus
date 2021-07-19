library(tidyverse)
library(vegan)
library(ggsci)

virus_list = read_csv("data/mammal_virus_list.csv")
host_list = read_csv("data/host_occurence.csv")
host_info = read_csv("data/host_metadata.csv")

host_merged = inner_join(
  host_info,
  host_list,
  by = "Mammal_species"
)

host_comm = host_merged %>%
  select(Jingmen, Longquan, Wenzhou, Wufeng) %>%
  t()
colnames(host_comm) = host_merged$Mammal_species

host_list_rarefied = host_comm %>%
  rrarefy(621) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Mammal_species = rownames(.)) %>%
  as_tibble() %>%
  left_join(
    host_info,
    by = "Mammal_species"
  )

## ############ PART1
all_hosts = host_info$Mammal_species
rodents = filter(host_info, Mammal=="rodents")$Mammal_species
bats = filter(host_info, Mammal=="bats")$Mammal_species
shrews = filter(host_info, Mammal=="shrews")$Mammal_species

calc_alpha_div = function(host, method, tit, nsample = 20) {
  sample_specnum = matrix(nrow=nsample, ncol=4)
  for (i in 1:nsample) {
    tmp = host_comm %>%
      rrarefy(620) %>%
      as_tibble() %>%
      select(all_of(host)) %>%
      method()
    sample_specnum[i,] = tmp
  }
  print(sample_specnum)
  colnames(sample_specnum) = c("Jingmen", "Longquan", "Wenzhou", "Wufeng")
  
  df_specnum = sample_specnum %>%
    t() %>%
    as.data.frame() %>%
    mutate(site_name = rownames(.)) %>%
    pivot_longer(cols = 1:nsample, names_to="nsample")
  
  ggplot(df_specnum) + 
    geom_bar(aes(x=site_name, y=value, fill=site_name), stat = "identity") +
    scale_fill_uchicago() +
    my_ggtheme + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
    ylab("Host Richness") + 
    xlab("") + 
    ggtitle(tit)
}

# Spec composition

r_comm = host_comm %>% rrarefy(621)
r_ = r_comm %>%
  as.data.frame() %>%
  select(all_of(rodents)) %>%
  rowSums()

b_ = r_comm %>%
  as.data.frame() %>%
  select(all_of(bats)) %>%
  rowSums()

s_ = r_comm %>%
  as.data.frame() %>%
  select(all_of(shrews)) %>%
  rowSums()

tmp = rbind(rodents = r_, bats = b_, shrews = s_) %>%
  as.data.frame() %>%
  mutate(Mammal = rownames(.)) %>%
  pivot_longer(cols=1:4, names_to="site_name")

ggplot(tmp, aes(x=site_name, y=value, fill=Mammal)) + 
  geom_bar(position="stack", stat = "identity")

### 
r_comm %>%
  as.data.frame() %>%
  select(all_of(rodents)) %>%
  vegdist(method = "jaccard", binary=TRUE)

r_comm %>%
  as.data.frame() %>%
  select(all_of(bats)) %>%
  vegdist(method = "jaccard", binary=TRUE)

r_comm %>%
  as.data.frame() %>%
  select(all_of(shrews)) %>%
  vegdist(method = "jaccard", binary=TRUE)

order(r_comm, decreasing = TRUE)

### END OF HOST DIVERSITY ANALYSIS

### VIRAL DIVERSITY

virus_table = virus_list %>%
  select(Viral_Family_Genus_Species, Total_Abundance, Host_Name) %>%
  pivot_wider(
    names_from = Viral_Family_Genus_Species, 
    values_from = Total_Abundance, 
    values_fill = 0)

wf = startsWith(virus_table$Host_Name, "Wufeng ")
jm = startsWith(virus_table$Host_Name, "Jingmen ")
lq = startsWith(virus_table$Host_Name, "Longquan ")
wz = startsWith(virus_table$Host_Name, "Wenzhou ")
site_names = rep("", nrow(virus_table))
site_names[wf] = "Wufeng"
site_names[jm] = "Jingmen"
site_names[lq] = "Longquan"
site_names[wz] = "Wenzhou"

virus_table$site_name = site_names

virus_table$Host_Name = virus_table$Host_Name %>%
  str_replace_all( "Wufeng ", "") %>%
  str_replace_all( "Jingmen ", "") %>%
  str_replace_all( "Longquan ", "") %>%
  str_replace_all( "Wenzhou ", "")

tr = virus_table %>%
  left_join(host_info, by=c("Host_Name"="Mammal_species")) %>%
  select(-Host_Name, -Mammal_genus) %>%
  group_by(site_name, Mammal) %>%
  summarise_all(function(vec){as.integer(sum(vec > 0)>=1)}) %>%
  ungroup() %>%
  pivot_longer(cols=3:ncol(.), names_to="virus_species") %>%
  select(-virus_species) %>%
  group_by(site_name, Mammal) %>%
  summarise_all(sum)

tmp = virus_table %>%
  left_join(host_info, by=c("Host_Name"="Mammal_species")) %>%
  select(-Host_Name, -Mammal_genus) %>% group_by(Mammal, site_name) %>% summarise(n=n())

inner_join(tr, tmp, by=c("site_name","Mammal")) %>%
  mutate(r = value/n)


##########
tmp = virus_table %>%
  left_join(host_info, by=c("Host_Name"="Mammal_species")) %>%
  select(-Mammal_genus) %>%
  group_by(site_name, Host_Name) %>%
  summarise_all(function(vec){as.integer(sum(vec > 0)>=1)}) %>%
  ungroup() %>%
  pivot_longer(cols=3:ncol(.), names_to="virus_species") %>%
  select(-virus_species) %>%
  group_by(site_name, Host_Name) %>%
  summarise_all(sum) %>%
  left_join(host_info, by=c("Host_Name"="Mammal_species"))

##########
tmp = virus_table %>%
  select(-site_name, -Host_Name) 

tmp = host_comm %>% 
  rrarefy(620) %>% 
  as.data.frame() %>%
  select(all_of(rodents))

hb = tmp[1,] + tmp[4,]
zj = tmp[2,] + tmp[3,]
rbind(hb, zj) %>%
  specnumber()

ha = tmp[2,] + tmp[4,]
la = tmp[1,] + tmp[3,]
rbind(ha, la) %>%
  specnumber()

for (i in 1:23) {
  title = str_glue(virus_table[i,]$Host_Name, " @ ", virus_table[i,]$site_name)
  
  k = sum(as.double(tmp[i,]) > 0)
  y = as.double(sort(as.data.frame(tmp[i,]), decreasing=T)[1:k])
  x = colnames(tmp)[order(as.double(tmp[i,]), decreasing=T)][1:k]
  ggplot() +
    geom_bar(aes(x=factor(x, levels=x),y=y*100000), stat="identity") +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
    ylab("RPM") + 
    xlab("") +
    ggtitle(title)
  ggsave(str_glue(title, ".pdf"))
}

