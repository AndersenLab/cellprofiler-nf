#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)
library(purrr)
library(data.table)

#==============================================================================#
# Arguments
#==============================================================================#
# 1 - full path to project directory
# 2 - the path to the well_mask - HARDCODE NOW
# 3 - the group argument from main.nf - default is "plate,well"
# 4 - the edited pipeline path
# 5 - the out path
args <- commandArgs(trailingOnly = TRUE)

#==============================================================================#
# Make Metadata NEEDS TO BE ADAPATBLE TO MULTIPLE WAVELENGTHS
#==============================================================================#
projDir <- args[1]
projName <- stringr::str_extract(projDir, pattern = "([^/]+$)")

raw_imagesDir <- paste0(projDir, "/raw_images")

# parse file names from directory - need wavelength in file name
meta1 <- tibble::tibble(file = list.files(path = raw_imagesDir),
                        file_path = list.files(path = raw_imagesDir, full.names = T)) %>%
  dplyr::filter(is.na(stringr::str_match(file, "(?:[0-9]+-[a-z|A-Z|0-9]+-p[0-9]+-m[0-9]+[x|X]_[A-Z][0-9]+_w[0-9])(_thumb)?(?:|[:word:]{36})\\.([a-z]{3}|[A-Z]{3})")[,2])) %>%
  dplyr::mutate(copy = stringr::str_match(file, "([0-9]+-[a-z|A-Z|0-9]+-p[0-9]+-m[0-9]+[x|X]_[A-Z][0-9]+_w[0-9])(?:|[:word:]{36})")[,2]) %>%
  dplyr::mutate(file = stringr::str_c(copy, ".TIF")) %>%
  #dplyr::mutate(copy = file) %>%
  tidyr::separate(col = copy, into = c("date","exp","plate","mag"), sep = "-") %>%
  tidyr::separate(col = mag, into = c("mag","well", "wave"), sep = "_") %>%
  #tidyr::separate(col = wave, into = c("wave","TIF"), sep = "[.]") %>%
  #dplyr::select(-TIF) %>%
  dplyr::mutate(row = stringr::str_extract(well, pattern = "[A-Z]"),
                col = stringr::str_extract(well, pattern = "[0-9][0-9]"),
                Image_PathName_wellmask_98.png = stringr::str_replace(args[2], pattern = "([^/]+$)", replacement = ""),
                Image_FileName_wellmask_98.png = stringr::str_extract(args[2], pattern = "([^/]+$)"))

# num of wavelengths - add logic for how to make metadata from multiple wavelengths
n_wave <- length(unique(meta1$wave))

# add group
groups <- stringr::str_split(args[3], pattern = ",")[[1]]
meta1$group <- apply( meta1[, groups], 1, paste, collapse = "_")
write.table(meta1, file = glue::glue("groups.tsv"), quote=FALSE, sep='\t', row.names = F)

# add image types and set metadata names - hardcode image names - needs to be flexible for multiple pipeline profiles
meta2 <- meta1 %>%
  tidyr::pivot_wider(names_from = wave, values_from = c(file, file_path)) %>%
  dplyr::rename(Image_FileName_RawBF = file_w1,
                Image_PathName_RawBF = file_path_w1,  
                Image_FileName_RawRFP = file_w2,
                Image_PathName_RawRFP = file_path_w2) %>%
  dplyr::mutate(Image_PathName_RawRFP = stringr::str_replace(Image_PathName_RawRFP, pattern = "([^/]+$)", replacement = ""),
                Image_PathName_RawBF = stringr::str_replace(Image_PathName_RawBF, pattern = "([^/]+$)", replacement = "")) %>%
  dplyr::select(Metadata_Experiment = exp,
                Metadata_Date = date,
                Metadata_Plate = plate,
                Metadata_Well = well,
                #Metadata_Column = col,
                #Metadata_Row = row,
                Metadata_Group = group,
                Metadata_Magnification = mag,
                Image_FileName_RawBF,
                Image_PathName_RawBF,
                Image_FileName_RawRFP,
                Image_PathName_RawRFP,
                Image_FileName_wellmask_98.png,
                Image_PathName_wellmask_98.png)

write.table(meta2, file = glue::glue("metadata.csv"), quote=FALSE, sep=',', row.names = F)

#==============================================================================#
# Make groups.tsv file for runCP
#==============================================================================#
gs <- meta2 %>%
  dplyr::distinct(Metadata_Group, .keep_all=T) %>%
  dplyr::mutate(group = paste0("Metadata_Group=", Metadata_Group),
                pipeline = args[4],
                output = paste0(args[5], "/CP_output/", Metadata_Group)) %>%
  dplyr::select(group:output)

write.table(gs, file = glue::glue("groups.tsv"), quote=FALSE, sep='\t', row.names = F)

#==============================================================================#
# Make dirs for CP output
#==============================================================================#
for(i in unique(gs$output)){
  dir.create(i, recursive = T)
}
