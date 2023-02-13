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
#args <- c("/projects/b1059/projects/Tim/cellprofiler-nf/debug/20220501_toxinDebug/Analysis-20230126")
args <- commandArgs(trailingOnly = TRUE)

#==============================================================================#
# Read CP output data
#==============================================================================#
# get the output for each model
dir <- glue::glue("{args[1]}/processed_data")

# get the wormobject data too
PrimaryObjects <- readr::read_csv(file = dir %>% fs::dir_ls(regexp = "WormObjects\\.csv$")) %>%
  dplyr::select(FileName_RawBF, ObjectNumber, AreaShape_Area:AreaShape_Solidity) %>%
  rename_all(.vars = 2:27, ~ paste0("po_", .x))

# read in files and manipulate with  
model_df_raw <- dir %>%
  fs::dir_ls(regexp = "_NonOverlappingWorms\\.csv$") %>% # find paths to csvs in dir
  purrr::map_dfr(readr::read_csv, .id = "model") %>%
  dplyr::mutate(Metadata_Date = as.integer(Metadata_Date), 
                model = stringr::str_remove_all(fs::path_file(as.character(model)), pattern = "_NonOverlappingWorms|.csv"),
                model = paste0(model, ".model.outputs")) %>%
  dplyr::arrange(model, ImageNumber)

# join the model outputs with the primary object shape data
model_df <- model_df_raw %>%
  dplyr::left_join(PrimaryObjects, by = c("FileName_RawBF" = "po_FileName_RawBF", "Parent_WormObjects" = "po_ObjectNumber"))

# split to list
model_df_list <- split.data.frame(model_df, model_df$model)

# export list items to global env
lapply(seq_along(model_df_list), function(i) assign(names(model_df_list)[i], model_df_list[[i]], envir = .GlobalEnv))

# save as R.data
proj_name <- stringr::str_extract(args[1], pattern = "[^\\/]+(?=(?:\\/[^\\/]+){1}$)")
run_stamp <- stringr::str_extract(args[1], pattern = "([^/]+$)")
save(list = c(ls(pattern = "model.outputs|wormobj")),
     file = glue::glue("{args[1]}/processed_data/{proj_name}_{run_stamp}.RData"))

# clean up extra CP_outputs for now
system(command = glue::glue("if [ -d {args[1]}/CP_output ]; then rm -Rf {args[1]}/CP_output; fi"))