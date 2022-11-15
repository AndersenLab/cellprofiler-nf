#!/usr/bin/env -S Rscript --vanilla
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(purrr)
library(data.table)

#==============================================================================#
# Arguments
#==============================================================================#
# 1 - A list of input files
# 2 - the path to the well_mask - HARDCODE NOW
# 3 - the group argument from main.nf - default is "plate,well"
# 4 - the edited pipeline path
# 5 - the out path
args <- commandArgs(trailingOnly = TRUE)

#==============================================================================#
# Make Metadata NEEDS TO BE ADAPATBLE TO MULTIPLE WAVELENGTHS
#==============================================================================#
# parse file names from directory - need wavelength in file name
meta1 <- read_delim(
  args[1], 
  col_names = FALSE, 
  delim = "\t") %>% 
  select(file_path = X1) %>% 
  extract(file_path, into = "file", remove = FALSE, regex = ".*/(.*)$") %>% 
  extract(file, 
          remove = FALSE, 
          regex = "^(.*)-(.*)-(.*)-(.*)_(.*)\\.(.*)$", 
          into = c("date","exp","plate","mag","well","TIF")) %>% 
  select(-TIF) %>%
  dplyr::mutate(row = stringr::str_extract(well, pattern = "[A-Z]"),
                col = stringr::str_extract(well, pattern = "[0-9][0-9]"),
                Image_PathName_wellmask_98.png = stringr::str_replace(args[2], pattern = "([^/]+$)", replacement = ""),
                Image_FileName_wellmask_98.png = stringr::str_extract(args[2], pattern = "([^/]+$)"))

# add group
groups <- stringr::str_split(args[3], pattern = ",")[[1]]
meta1$group <- apply( meta1[, groups], 1, paste, collapse = "_")

# add image types and set metadata names - hardcode image names - needs to be flexible for multiple pipeline profiles
meta2 <- meta1 %>%
  dplyr::mutate(Image_PathName_RawBF = stringr::str_replace(file_path, pattern = "([^/]+$)", replacement = "")) %>%
  dplyr::select(Metadata_Experiment = exp,
                Metadata_Date = date,
                Metadata_Plate = plate,
                Metadata_Well = well,
                Metadata_Group = group,
                Metadata_Magnification = mag,
                Image_FileName_RawBF = file,
                Image_PathName_RawBF,
                Image_FileName_wellmask_98.png,
                Image_PathName_wellmask_98.png)

write.table(meta2, file = "metadata.csv", quote=FALSE, sep=',', row.names = F)