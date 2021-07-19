library(tidyverse)

path_to_virus_list = function() {
  "data/mammal_virus_list.csv"
}

read_virus_list = function(filename) {
  read_csv(filename)
}

# 病毒器官嗜性CA图
create_plot_ca_organophilism = function(virus_list) {
  virus_list %>%
    select(c("lung_Abundance",
             "liver_Abundance",
             "spleen_Abundance",
             "kindey_Abundance",
             "feces_Abundance")) %>%
    vegan::cca() %>%
    plot()
}

"Viral_Family"
"Viral_Genus"
"Mammal"
"Host_Name"

file_virus_list = path_to_virus_list()
tb_virus_list = read_virus_list(file_virus_list)

ca = tb_virus_list %>%
  select(c("lung_Abundance",
           "liver_Abundance",
           "spleen_Abundance",
           "kindey_Abundance",
           "feces_Abundance")) %>%
  vegan::cca()

ggplot() + 
  geom_point(aes(x=ca$CA$u[,1], y=ca$CA$u[,2], color=tb_virus_list$Viral_Family)) + 
  geom_segment(aes(x=0,y=0,xend=ca$CA$v[,1], yend=ca$CA$v[,2],),arrow=arrow()) + 
  geom_text(aes(x=ca$CA$v[,1], y=ca$CA$v[,2], label=rownames(ca$CA$v)))
