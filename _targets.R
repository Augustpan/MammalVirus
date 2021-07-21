library(targets)
source("R/functions.R")

#tar_option_set(debug = "xxx")

list(
  tar_target(file_virus_list, getPathToVirusList(), format = "file"),
  tar_target(file_host_list, getPathToHostList(), format = "file"),
  tar_target(file_host_info, getPathToHostInfo(), format = "file"),
  
  tar_target(tb_virus_list, read_csv(file_virus_list)),
  tar_target(tb_host_list, read_csv(file_host_list)),
  tar_target(tb_host_info, read_csv(file_host_info)),
  
  tar_target(tb_host_comm, buildHostCommunityTable(tb_host_list, tb_host_info)),
  tar_target(tb_specnum_all_hosts, calcHostAlphaDiversity(tb_host_comm, tb_host_info, "all", vegan::specnumber)),
  tar_target(gg_specnum_all_hosts, createPlotHostAlphaDiversity(tb_specnum_all_hosts))
)

#tar_make(callr_function = TRUE)