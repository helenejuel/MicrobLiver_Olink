# Load packages
source(here::here("R/package-loading.R"))

# Import data
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
renamed_olink <- olink %>%
    rename(MCP3 = `MCP-3`,
           LAP_TGFb = `LAP TGF-beta-1`,
           IL17C = `IL-17C`,
           MCP1 = `MCP-1`,
           IL17A = `IL-17A`,
           IL20RA = `IL-20RA`,
           IL2RB = `IL-2RB`,
           IL1a = `IL-1 alpha`,
           TGFa = `TGF-alpha`,
           MCP4 = `MCP-4`,
           FGF23 = `FGF-23`,
           IL10RA = `IL-10RA`,
           FGF5 = `FGF-5`,
           MMP1 = `MMP-1`)
# TODO: include all colnames with special characters

# TODO: set LOD and change values <LOD

# TODO: change time_point and sample_type to factor


# filter cohorts
alco <- renamed_olink %>%
    filter(cohort == "ALCO")
ald <- renamed_olink %>%
    filter(cohort == "ALD")

#



# Plotting ----------------------------------------------------------------

theme_set(theme_bw())

# histogram
alco %>%
    ggplot(aes(x = IL1a)) +
    geom_density()
