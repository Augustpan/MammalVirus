library(targets)
source("R/functions.R")

tar_option_set(debug = "xxx")

list(
  tar_target(
    file_virus_list,
    path_to_virus_list(),
    format = "file"
  ),
  tar_target(
    tb_virus_list,
    read_virus_list(file_virus_list)
  ),
  tar_target(
    
  )
)

tar_make(callr_function = TRUE)