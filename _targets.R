library(targets)
source("R/functions.R")

#tar_option_set(debug="tb_virus_list_clean")
list(
  # manage path to raw data
  tar_target(file_virus_list, getPathVirusList(), format = "file"),
  tar_target(file_host_list, getPathHostList(), format = "file"),
  tar_target(file_host_info, getPathHostInfo(), format = "file"),
  tar_target(file_host_spe_tree, getPathHostSpeciesTree(), format = "file"),
  tar_target(file_host_ord_tree, getPathHostOrderTree(), format = "file"),
  tar_target(file_virus_segment, getPathVirusSegment(), format = "file"),
  tar_target(file_virus_fam_pre_mat, getPathVirusFamPreMat(), format = "file"),
  tar_target(file_host_ord_pre_mat, getPathHostOrdPreMat(), format = "file"),
  tar_target(file_seq_stat, getPathSeqStat(), format = "file"),
  
  # read raw data
  tar_target(tb_virus_list, read_csv(file_virus_list, show_col_types = FALSE)),
  tar_target(tb_host_list, read_csv(file_host_list, show_col_types = FALSE)),
  tar_target(tb_host_info, read_csv(file_host_info, show_col_types = FALSE)),
  tar_target(tb_virus_segment, read_tsv(file_virus_segment, show_col_types = FALSE)),
  tar_target(tree_host_spe, read.nexus(file_host_spe_tree)),
  tar_target(tree_host_ord, read.tree(file_host_ord_tree)),
  tar_target(mat_virus_fam_pre_mat, readRDS(file_virus_fam_pre_mat)),
  tar_target(mat_host_ord_pre_mat, readRDS(file_host_ord_pre_mat)),
  tar_target(tb_seq_stat, read_csv(file_seq_stat , show_col_types = FALSE)),
  
  # data cleaning and transformation
  tar_target(tb_virus_list_clean, cleanVirusList(tb_virus_list, tb_virus_segment)),

  # analysis
  tar_target(tb_specnum_hosts, calcHostAlphaDiversity(tb_host_list, tb_host_info)),
  tar_target(tb_virus_tab_sorted, sortVirusTableCA(tb_virus_list)),
  tar_target(list_virus_annotation, summariseVirusAnnotation(tb_virus_list)),
  tar_target(tb_full_model_list, makeFullModelList(tb_virus_list_clean, mat_virus_fam_pre_mat, mat_host_ord_pre_mat)),
  #tar_target(tb_fitted_model, fitModel(tb_full_model_list, tb_virus_list_clean), pattern = map(tb_full_model_list)),
  #tar_target(tb_permanova_virus_op, performPERMANOVAVirusOp(tb_virus_list)),
  
  # create plots
  tar_target(gg_shared_hosts, createPlotSharedHostSpecies(tb_host_list, tb_host_info)),
  tar_target(gg_virus_op_ca, createPlotCAVirusOrganophilism(tb_virus_list))
  
)

