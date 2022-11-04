#!/usr/bin/env -S Rscript --vanilla
library(fs)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(purrr)

#==============================================================================#
# Arguments
#==============================================================================#
# 1 - out directory path
# 2 - project name
# 3 - run stamp
args <- commandArgs(trailingOnly = TRUE)

#==============================================================================#
# Read CP output data
#==============================================================================#
# get the output for each model
dir <- "processed_data"

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
proj_name <- args[1]
run_stamp <- args[2]
save(list = c(ls(pattern = "model.outputs")),
     file = paste0("processed_data/", args[1], "_", args[2], ".RData"))