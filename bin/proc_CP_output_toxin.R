#!/usr/bin/env Rscript
library(fs)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)
library(purrr)

#==============================================================================#
# Arguments
#==============================================================================#
# 1 - out directory path
args <- commandArgs(trailingOnly = TRUE)

#==============================================================================#
# Read CP output data
#==============================================================================#
# get the output for each model
dir <- glue::glue("{args[1]}/processed_data")

# read in files and manipulate with  
model_df <- dir %>%
  fs::dir_ls(regexp = "\\.csv$") %>% # find paths to csvs in dir
  purrr::map_dfr(readr::read_csv, .id = "model") %>%
  dplyr::mutate(Metadata_Date = as.integer(Metadata_Date), 
                model = stringr::str_remove(fs::path_file(as.character(model)), pattern = ".csv"),
                model = paste0(model, ".model.outputs")) %>%
  dplyr::arrange(model, ImageNumber)

# split to list
model_df_list <- split.data.frame(model_df, model_df$model)

# export list items to global env
lapply(seq_along(model_df_list), function(i) assign(names(model_df_list)[i], model_df_list[[i]], envir = .GlobalEnv))

# save as R.data
proj_name <- stringr::str_extract(args[1], pattern = "[^\\/]+(?=(?:\\/[^\\/]+){1}$)")
run_stamp <- stringr::str_extract(args[1], pattern = "([^/]+$)")
save(list = c(ls(pattern = "model.outputs")),
     file = glue::glue("{args[1]}/processed_data/{proj_name}_{run_stamp}.RData"))

# clean up extra CP_outputs for now
system(command = glue::glue("if [ -d {args[1]}/CP_output ]; then rm -Rf {args[1]}/CP_output; fi"))