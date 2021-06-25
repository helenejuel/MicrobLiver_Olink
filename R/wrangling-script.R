# Load packages
source(here::here("R/package-loading.R"))

# Import data
olink <- readxl::read_excel(here::here("data-raw/20202249_Juel_NPX.xlsx"), sheet = "results")
meta <- readxl::read_excel(here::here("data-raw/20202249_Juel_NPX.xlsx"), sheet = "assay_meta")

# Save dfs in data folder
usethis::use_data(olink, overwrite = T)
usethis::use_data(meta, overwrite = T)

