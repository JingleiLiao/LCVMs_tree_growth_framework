# =====================================================
# 01_data_preparation.R
#
# Data preparation for hierarchical climate-growth analysis
#
# Prepare:
# 1. Tree-ring chronology
# 2. Regional climate variables
# 3. Large-scale climate variability modes
#
# =====================================================


# =====================================================
# 1. Load packages
# =====================================================

library(dplyr)
library(tidyr)
library(tibble)


# =====================================================
# 2. Define data paths
# =====================================================
#
# Users should modify these paths according to
# their local data organization.
#
# =====================================================


data_path <- "your_data_directory"


tree_file <- file.path(
  data_path,
  "tree_ring_chronology.csv"
)


climate_file <- file.path(
  data_path,
  "regional_climate.csv"
)


lcvm_file <- file.path(
  data_path,
  "large_scale_climate_modes.csv"
)



# =====================================================
# 3. Load datasets
# =====================================================


# Tree-ring chronology dataset
#
# Required columns:
# year : calendar year
# lon  : longitude
# lat  : latitude
# std  : standardized tree-ring chronology


STD <- read.csv(
  tree_file,
  stringsAsFactors = FALSE
)



# Regional climate dataset
#
# Required columns:
# year
# lon
# lat
# PG   : precipitation
# TG   : temperature
# SPEIG: drought index


T <- read.csv(
  climate_file,
  stringsAsFactors = FALSE
)



# Large-scale climate variability modes
#
# Required variables:
# ENSO
# PDO
# NAO
# AO
# Sea ice


AI <- read.csv(
  lcvm_file,
  stringsAsFactors = FALSE
)



# =====================================================
# 4. Basic data checking
# =====================================================


# Check data structure

str(STD)

str(T)

str(AI)



# Check missing values

colSums(is.na(STD))

colSums(is.na(T))

colSums(is.na(AI))