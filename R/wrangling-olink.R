# Load packages
source(here::here("R/package-loading.R"))


# Load raw data -> olink4 -------------------------------------------------

# Load data from raw NPX data
olink <- OlinkAnalyze::read_NPX(here::here("data-raw/20202249_Juel_NPX.xlsx"))

# Let the wrangling begin!
str(olink)

# Fix issue with TIPS naming
# replace typo TIEP80EC	--> TIPS80EC
olink2 <- olink %>%
    mutate(SampleID = sub("TIEP80EC", "TIPS80EC", SampleID))
# replace typo TIPSKC/93 --> TIPS93
olink2 <- olink2 %>%
    mutate(SampleID = sub("TIPSKC/93", "TIPS93", SampleID))

# Add columns for cohort and pt_ID to be able to extract each cohort
# I have a sampleID column with values in the format "cohort-patientID-timepoint-sampletype" that I ultimately need to break down to 4 separate columns.
# The cohort is 2-7 UPPERCASE letters, which is written like this "^[A-Z]{2,7}"
# No separator, then patientID is 1-4 numbers, written like this: "[0-9]{1,4}"
# Only some samples have timepoint and sampletype, and they have different patterns depending on the cohort, so I think the easiest will be to first extract cohort and patientID, and then extract the last 2 for each cohort separately.
olink2 <- olink2 %>%
    mutate(cohort = sub("[0-9]{1,4}.*$", "", SampleID)) %>% # substitute 1-4 numbers and everything after this (.*$) with nothing
    mutate(cohort = sub("_.*$", "", cohort)) # substitute anything after a _ with nothing to remove long control names and just keep "CONTROL"

# Keep in mind that * means "match at least 0 times" while + means "match at least 1 time". For some reason the code using .+$ did not replace the single number in the ppts with a single cifre ID.

# Transform to factor to check correct number of cohorts
olink3 <- olink2
olink3$cohort <- as.factor(olink3$cohort)
summary(olink3$cohort) # if transformed to factor for check: 9 different cohorts, including "Bridge" and "CONTROL"

# change cytokine names
# first remove all "-"
# then change spaces to "_"
# then rename alpha/beta/gamma to a/b/g
olink3 <- olink2 %>%
    mutate(Assay = gsub("-", "", Assay)) %>%
    mutate(Assay = gsub(" ", "_", Assay)) %>%
    mutate(Assay = gsub("alpha", "a", Assay)) %>%
    mutate(Assay = gsub("beta", "b", Assay)) %>%
    mutate(Assay = gsub("gamma", "g", Assay))

# count (remove) samples with QC warning
olink3$QC_Warning <- as.factor(olink3$QC_Warning)
# str(olink3)
summary(olink3$QC_Warning) #2576 have warning

# To remove all samples with QC warning (will be done for each cohort separately to make it easier to study the warning samples)
# olink3 <- olink3 %>%
#     filter(QC_Warning == "Pass")

# change values < LOD of that assay to 50% LOD
olink4 <- olink3 %>%
    mutate(corr_NPX = if_else(LOD > NPX, LOD/2, NPX))

# save dataset in data folder
usethis::use_data(olink4, overwrite = T)

# Load dataset
load(here::here("data/olink4.rda"))

# Olink ID list
cytokines <- levels(as.factor(olink4$Assay)) # protein abbreviations
cytokine_IDs <- levels(as.factor(olink4$OlinkID))


# Gastric Bypass ----------------------------------------------------------
# Create columns for pt_ID, time_point, sample_type
MLGB <- olink4 %>%
  filter(cohort == "MLGB") %>%
  mutate(pt_ID = sub("^[A-Z]{4}", "", sub("-.+$", "", SampleID))) %>%
  mutate(time_point = sub("^.+-", "", SampleID)) %>%
  mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID))

# Remove duplicate sample
str(as.factor(MLGB$sample_type)) # we can see that the duplicate sample has sample_type "M-d"
# MLGB %>%
#     filter(sample_type == "M-d") #92 duplicates (1 sample)
MLGB <- MLGB %>%
  filter(sample_type != "M-d") # goes from 14812 to 14720 obs

# remove samples with QC warning
# MLGB %>%
#     filter(QC_Warning != "Pass") # 92 obs
MLGB <- MLGB %>%
  filter(QC_Warning == "Pass") # goes from 14720 to 14628 obs

# Remove sample_type column
MLGB$sample_type <- NULL

# Change time_point and pt_ID to factor
MLGB <- MLGB %>%
  mutate(across(c(time_point, pt_ID), factor))

# Change order of timepoints
MLGB$time_point <- ordered(MLGB$time_point, levels = c("BL", "3M", "12M"))

# str(MLGB)

# save dataset in data folder
usethis::use_data(MLGB, overwrite = T)

# Load dataset
# load(here::here("data/MLGB.rda"))

# Alcochallenge -----------------------------------------------------------

# Create columns for pt_ID, time_point, sample_type
ALCO <- olink4 %>%
    filter(cohort == "ALCO") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", sub("-.+$", "", SampleID))) %>%  # substitute 2-7 letters and anything after a "-" with nothing
    mutate(time_point = sub("^.+-T", "", sub("S.+$", "", SampleID))) %>% # add timepoint
    mutate(sample_type = sub("^.+S", "", SampleID)) # add sample_type (peripheral or liver vein)

# Remove duplicate sample
# ALCO %>%
#     filter(sample_type == "1-d") #92 duplicates (1 sample)
ALCO <- ALCO %>%
    filter(sample_type != "1-d") # goes from 22172 to 22080 obs

# remove samples with QC warning
# ALCO_fail <- ALCO %>%
#     filter(QC_Warning != "Pass") # 552 obs - all pt 106
# ALCO_fail %>%
#     filter(pt_ID != "106")
ALCO <- ALCO %>%
    filter(QC_Warning == "Pass") # goes from 22080 to 21528 obs


# Add subgroup column: 100=ALD, 200=NAFLD, 300=healthy
ALCO$pt_ID <- as.numeric(ALCO$pt_ID)
ALCO <- ALCO %>%
    mutate(subgroup = if_else(pt_ID > 300, "Healthy", if_else(pt_ID > 200, "NAFLD", "ALD")))

# Change pt_ID, time_point, sample_type and subgroup to factor
ALCO$pt_ID <- as.factor(ALCO$pt_ID) # 39 indv - 106 removed in QC warning
ALCO$time_point <- as.factor(ALCO$time_point) # 3 TPs
ALCO$sample_type <- as.factor(ALCO$sample_type) # 2 sample types
ALCO$subgroup <- as.factor(ALCO$subgroup) # 3 subgroups

# Change order of subgroups
ALCO$subgroup <- ordered(ALCO$subgroup, levels = c("Healthy", "NAFLD", "ALD"))

str(ALCO)

# save dataset in data folder
usethis::use_data(ALCO, overwrite = T) # not done yet

# Load dataset
# load(here::here("data/ALCO.rda"))



# Rifaximin ---------------------------------------------------------------

RFX <- olink4 %>%
    filter(cohort == "RFX") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", sub("-.+$", "", SampleID))) %>%
    mutate(time_point = sub("^.+-V", "", SampleID))

# Remove duplicate sample
RFX <- RFX %>%
    filter(time_point != "1-d")

# Change pt_ID, time_point and sample_type to factor
RFX$pt_ID <- as.factor(RFX$pt_ID) # 136 indv
RFX$time_point <- as.factor(RFX$time_point) # 3 TPs
str(RFX)

RFX %>% filter(time_point == 1) # 12,328 rows = 134 ppts
RFX %>% filter(time_point == 2) # 11,500 rows = 125 ppts
RFX %>% filter(time_point == 3) # 6,532 rows = 71 ppts



# TIPS --------------------------------------------------------------------

TIPS <- olink4 %>%
    filter(cohort == "TIPS") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", SampleID)) %>%
    mutate(pt_ID = sub("[A-Z]{2}", "", pt_ID)) %>% #remove the 2 letters
    mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID))




# GALA --------------------------------------------------------------------

### Load olink dataset
load(here::here("data/olink4.rda"))

GALA_olink <- olink4 %>%
    filter(cohort == "ALD" | cohort == "HP") %>% #50692 obs = 551 samples
    mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID)) # duplicate samples are marked as "-d"
# str(as.factor(GALA_olink$sample_type)) # check we have the duplicate samples as "-d"

# Remove duplicate samples
GALA_olink <- GALA_olink %>%
    filter(sample_type != "-d")  # goes from 50692 to 50416 obs = 276 duplicates (3 samples)

GALA_olink$sample_type <- NULL # remove sample_type column

# Assess samples with QC warning
# fail <- GALA_olink %>%
#     filter(QC_Warning != "Pass") # 1012 obs (12 samples)
# unique(fail$SampleID)

# remove samples with QC warning
GALA_olink <- GALA_olink %>%
    filter(QC_Warning == "Pass") # goes from 50416 to 49404 obs

# Remove extra 0's in SampleID HP004-HP096
GALA_olink <- GALA_olink %>%
  mutate(SampleID = sub("^HP[0]{1,2}", "HP", SampleID))
# print(unique(sort(GALA_olink$SampleID)))

# Remove cohort column, so it is not duplicate after merging
GALA_olink$cohort <- NULL



### Load phenotype data
# TODO: Change to new QC'ed phenotype dataset
GALA_pheno2 <- read_excel(here::here("data-raw/ALD_HP_outcome_HBJ.xlsx"))

# Change data columns to factor
GALA_pheno2 <- GALA_pheno2 %>%
  mutate(across(c(cohort, gender, fibrosis, inflam, abstinent, overuse, alcoholyears, abstinenceyears, liverrelated_event, hospInfection, allMortality, excess_drink_followup), factor)) # consider including kleiner, starts_with("nas")

# Change order of factors, if needed
GALA_pheno2$cohort <- ordered(GALA_pheno2$cohort, levels = c("HP", "ALD"))


# Change data columns to numeric
GALA_pheno2 <- GALA_pheno2 %>%
  mutate(across(c(cpa, meld, elf, te, packyears, height, weight, bmi, waist, hip, whr, hr, map, sbp, dbp, alt, ast, alk, bili, chol, crp, ggt, glc, hba1c, hdl, iga, igg, igm, ldl, leu, trigly, insulin, homair, cpeptid, proc3, days_to_LRE, days_to_hospInf, days_to_mort), as.numeric))

# Change column name from CBMR_ID to SampleID to allow merging with olink data
GALA_pheno2 <- GALA_pheno2 %>%
  rename(SampleID = CBMR_ID)

# str(GALA_pheno2)


### Match and merge Olink with phenotype data
GALA <- merge(GALA_olink, GALA_pheno2, by = "SampleID") # goes from 49404 to 49312 obs from GALA_olink = 1 ppt!
# Find missing sample
# print(setdiff(unique(GALA_olink$SampleID), GALA_pheno2$SampleID)) # ALD1302
# TODO: include phenotypes for ALD1302

# Add column for fibrosis level high or low
GALA <- GALA %>%
  mutate(te_fibrosis = if_else(te < 6, "low", "high"))
GALA$te_fibrosis <- as.factor(GALA$te_fibrosis)

# save dataset in data folder
usethis::use_data(GALA, overwrite = T)


### If using Nanna's phenotype dataset instead
# GALA_SIP_pheno <- read_excel(here::here("data-raw/ALD_HP_SIP_phenotypes.xlsx"), sheet = "small") #39 variables
# # str(GALA_SIP_pheno[,c(1:10)]) # Original phenotype columns converted to class chr
#
# # Remove PRS columns
# GALA_SIP_pheno <- GALA_SIP_pheno[,c(1:28)] #28 variables
#
# # Change gender, abstinent to factor
# # GALA_SIP_pheno$gender <- as.factor(GALA_SIP_pheno$gender) # old method single
#
# # old method multiple
# # for(i in c(2,4)) {
# #     GALA_SIP_pheno[, i] = as.factor(GALA_SIP_pheno[, i])
# #     }
#
# # With dplyr - note: just 'factor', not 'as.factor' (though both seem to work)
# GALA_SIP_pheno <- GALA_SIP_pheno %>%
#   mutate(across(c(gender, abstinent, kleinerFscore), factor))
#
# #Change age GALA_SIP_pheno[,3] and clinical parameters GALA_SIP_pheno[,c(5:23, 26, 28)] to numeric
# GALA_SIP_pheno <- GALA_SIP_pheno %>%
#   mutate(across(all_of(names(GALA_SIP_pheno)[c(3,5:26,28)]), as.numeric))
# # str(GALA_SIP_pheno)
#
# # filter for ALD and HP
# GALA_pheno <- GALA_SIP_pheno %>%
#     mutate(cohort = sub("[0-9]{1,5}", "", ID)) %>%  # substitute 1-5 numbers with nothing
#     filter(cohort == "ALD" | cohort == "HP") %>%
#     mutate(SampleID = ID) %>% # add SampleID column to be able to merge with Olink data
#     mutate(pt_ID = sub("[A-Z]{2,3}", "", ID)) #remove the 2-3 letters in a new pt_ID column to merge with SIP
# # str(as.factor(GALA_pheno$cohort)) # check: 2 levels left
#
# # Filter SIP samples to be able to include the ppts that are included in GALAXY
# SIP_pheno <- GALA_SIP_pheno %>%
#     mutate(cohort = sub("[0-9]{1,5}", "", ID)) %>%
#     filter(cohort == "SIP") %>%
#     mutate(pt_ID = sub("[A-Z]{3}", "", ID)) %>% #remove the 3 letters
#     mutate(SampleID = paste0("ALD", pt_ID))
#
# # rbind GALA_pheno with SIP
# GALA_SIP_pheno2 <- rbind(GALA_pheno, SIP_pheno) # back to 1187 individuals
# GALA_SIP_pheno2$cohort <- NULL # remove cohort column before merging, to prevent duplicate
#
# # match and merge Olink with phenotype data
# GALA <- merge(GALA_olink, GALA_SIP_pheno2, by = "SampleID") # goes from 49404 to 46276 obs from GALA_olink = 34 ppts
#
# # Figuring out which IDs are being dropped:
# GALA_olink_IDs <- sort(unique(GALA_olink$SampleID)) #537 samples
# # GALA_pheno_IDs <- sort(unique(GALA_pheno$SampleID)) #458 samples --> SIPHON-ALD samples removed from phenotype data! Fixed by matching ID number between GALA_olink and GALA_pheno for SIP only, but still some phenotypes missing or not matched
# # GALA_SIP_pheno2_IDs <- sort(unique(GALA_SIP_pheno2$SampleID)) #1181 unique SampleIDs of 1187 samples...?
# GALA_IDs <- sort(unique(GALA$SampleID)) #503 samples i.e. missing 34 samples from olink
#
# # Extract all the differing Sample IDs in the two vectors
# print(setdiff(GALA_olink_IDs, GALA_IDs))
# # These 34 are missing from Nanna's phenotype data, maybe because they were not genotyped or failed QC?
# Probably would have been easier to filter for ALD and HP separately, then match with phenotype data, and then rbind.


### Add outcome data if using dataset without outcome data
# # Load outcome data
# # TODO: Change to new QC'ed phenotype dataset
# GALA_outcome <- read_excel(here::here("data-raw/GALA-ALD_outcome_data.xlsx"))
# # str(GALA_outcome) # columns types messed up
#
# # Change id to integer
# GALA_outcome$id <- as.integer(GALA_outcome$id)
# # Change starts_with(days_to) columns to numeric
# GALA_outcome <- GALA_outcome %>%
#   mutate(across(starts_with("days_to"), as.numeric))
# # Change liver-related_event, hospInfection, allMortality, excess_drink_followup columns to factor
# GALA_outcome <- GALA_outcome %>%
#   mutate(across(c(liverrelated_event, hospInfection, allMortality, excess_drink_followup), factor))
# str(GALA_outcome) # check
#
# # Merge GALA with outcome data
# # Create SampleID column in outcome data
# GALA_outcome <- GALA_outcome %>%
#   mutate(SampleID = CBMR_ID)
# # Merge
# GALA_FU <- merge(GALA, GALA_outcome, by = "SampleID") # 34500 obs = 375 individuals -> dropped from 462 in outcome data and 537 in olink data?
# str(GALA_FU)
# GALA_FU_IDs <- sort(unique(GALA_FU$SampleID)) #375 samples
# # Compare the character vectors
# GALA_IDs %in% GALA_FU_IDs
# # TODO: extract the FALSE Sample IDs from the comparison
# GALA %>%
#   filter(SampleID == "ALD2693")
#
# # save dataset in data folder
# usethis::use_data(GALA_FU, overwrite = T)
