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
# args <- c("/projects/b1059/projects/Tim/cellprofiler-nf/projects/20220128_GWA09", "/projects/b1059/projects/Tim/cellprofiler-nf/input_data/well_masks/wellmask_98.png",
#          "plate,well", "/projects/b1059/projects/Tim/cellprofiler-nf/projects/20220128_GWA09/pipelines/pipeline.cppipe", "/projects/b1059/projects/Tim/cellprofiler-nf/projects/20220128_GWA09/CP_output")
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
  dplyr::filter(grepl("([0-9]+-[a-z|A-Z|0-9|_]+-p[0-9]+-m[0-9]+[X|x]_[A-Z][0-9]{2})(\\.TIF)$", file)) %>%
  dplyr::mutate(copy = file) %>%
  tidyr::separate(col = copy, into = c("date","exp","plate","mag"), sep = "-") %>%
  tidyr::separate(col = mag, into = c("mag","well"), sep = "_") %>%
  tidyr::separate(col = well, into = c("well","TIF"), sep = "[.]") %>%
  dplyr::select(-TIF) %>%
  dplyr::mutate(row = stringr::str_extract(well, pattern = "[A-Z]"),
                col = stringr::str_extract(well, pattern = "[0-9][0-9]"),
                Image_PathName_wellmask_98.png = stringr::str_replace(args[2], pattern = "([^/]+$)", replacement = ""),
                Image_FileName_wellmask_98.png = stringr::str_extract(args[2], pattern = "([^/]+$)"))

# add group
groups <- stringr::str_split(args[3], pattern = ",")[[1]]
meta1$group <- apply( meta1[, groups], 1, paste, collapse = "_")

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
