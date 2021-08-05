# Load packages
source(here::here("R/package-loading.R"))

# Import data
### OLD METHOD!!!###
olink <- readxl::read_excel(here::here("data-raw/20202249_Juel_NPX.xlsx"), sheet = "results")
meta <- readxl::read_excel(here::here("data-raw/20202249_Juel_NPX.xlsx"), sheet = "assay_meta")

# Save dfs in data folder
usethis::use_data(olink, overwrite = T)
usethis::use_data(meta, overwrite = T)

# Load data
load(here::here("data/olink.rda"))
load(here::here("data/meta.rda"))

# Let the wrangling begin!
str(olink)

# change cytokine names
# first remove all -
# renamed_olink <- rename_with(olink, ~ gsub("-", "", .x))
# then change spaces to _
# then rename alpha/beta/gamma to a/b/g
renamed_olink <- olink %>%
    rename_with(~ gsub("-", "", .x)) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(~ gsub("alpha", "a", .x)) %>%
    rename_with(~ gsub("beta", "b", .x)) %>%
    rename_with(~ gsub("gamma", "g", .x))
str(renamed_olink)
# code to change an individual colname:
# renamed_olink <- renamed_olink %>%
    # rename(LAP_TGFb = `LAP TGF-beta-1`)

# now rename cytokines in the metadata as well
renamed_meta <- meta %>%
    mutate(Assay = gsub("-", "", Assay)) %>%
    mutate(Assay = gsub(" ", "_", Assay)) %>%
    mutate(Assay = gsub("alpha", "a", Assay)) %>%
    mutate(Assay = gsub("beta", "b", Assay)) %>%
    mutate(Assay = gsub("gamma", "g", Assay))

# And rename column
renamed_meta <- renamed_meta %>%
    rename(missing_data_freq = `Missing Data freq.`)

