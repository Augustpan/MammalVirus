library(targets)
source("R/functions.R")

tar_option_set(debug="tb_permanova_virus_op")
list(
  # manage path to raw data
  tar_target(file_virus_list, getPathToVirusList(), format = "file"),
  tar_target(file_host_list, getPathToHostList(), format = "file"),
  tar_target(file_host_info, getPathToHostInfo(), format = "file"),
  tar_target(file_virus_list_annotated, getPathToAnnotatedVirusList(), format = "file"),
  
  # read raw data
  tar_target(tb_virus_list, read_csv(file_virus_list,show_col_types = FALSE)),
  tar_target(tb_host_list, read_csv(file_host_list,show_col_types = FALSE)),
  tar_target(tb_host_info, read_csv(file_host_info,show_col_types = FALSE)),
  tar_target(tb_virus_list_annotated, read_csv(file_virus_list_annotated,show_col_types = FALSE)),

  # analysis
  tar_target(tb_specnum_hosts, calcHostAlphaDiversity(tb_host_list, tb_host_info)),
  tar_target(tb_virus_tab_sorted, sortVirusTableCA(tb_virus_list)),
  tar_target(list_virus_annotation, summariseVirusAnnotation(tb_virus_list_annotated)),
  tar_target(tb_permanova_virus_op, performPERMANOVAVirusOp(tb_virus_list)),
  
  # create plots
  tar_target(gg_shared_hosts, createPlotSharedHostSpecies(tb_host_list, tb_host_info)),
  tar_target(gg_virus_op_ca, createPlotCAVirusOrganophilism(tb_virus_list))
  
)
