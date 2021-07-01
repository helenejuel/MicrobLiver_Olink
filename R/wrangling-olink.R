# Load packages
source(here::here("R/package-loading.R"))

# Alternative data loading from raw NPX data
olink2 <- OlinkAnalyze::read_NPX(here::here("data-raw/20202249_Juel_NPX.xlsx"))

# Let the wrangling begin!
str(olink2)

# Fix issue with TIPS naming
# replace typo TIEP80EC	--> TIPS80EC
olink3 <- olink2 %>%
    mutate(SampleID = sub("TIEP80EC", "TIPS80EC", SampleID))
# replace typo TIPSKC/93 --> TIPS93
olink3 <- olink3 %>%
    mutate(SampleID = sub("TIPSKC/93", "TIPS93", SampleID))

# Add columns for cohort and pt_ID to be able to extract each cohort
# I have a sampleID column with values in the format "cohort-patientID-timepoint-sampletype" that I ultimately need to break down to 4 separate columns.
# The cohort is 2-7 UPPERCASE letters, which I believe is written like this "^[A-Z]{2,7}"
# No separator, then patientID is 1-4 numbers, written like this: "[1-9]{1-4}"
# Only some samples have timepoint and sampletype, and they have different patterns depending on the cohort, so I think the easiest will be to first extract cohort and patientID, and then extract the last 2 for each cohort separately.
#grep selects all rows that have this pattern in the SampleID
#grep("^[A-Z]{2,7}", olink2$SampleID) # from the beginning, 2-7 uppercase letters for cohort
#grep("[1-9]{1,4}", olink2$SampleID) # 1-4 numbers for pt_ID
olink4 <- olink3 %>%
    mutate(cohort = sub("[0-9]{1,4}.*$", "", SampleID)) %>% # substitute 1-4 numbers and everything after this (.*$) with nothing
    mutate(cohort = sub("_.*$", "", cohort)) # substitute anything after a _ with nothing to remove long control names and just keep "CONTROL"

#Keep in mind that * means "match at least 0 times" while + means "match at least 1 time". For some reason the code using .+$ did not replace the single number in the ppts with a single cifre ID.

# olink4$cohort <- as.factor(olink4$cohort)

str(olink4)
summary(olink4$cohort) # if transformed to factor for check: 9 different cohorts, including "Bridge" and "CONTROL"

# TODO: Mark duplicate samples in a new column
# olink5 <- olink4 %>%
#     mutate(duplicate = grep("d.+$", SampleID)) #not working

# change cytokine names
# first remove all -
# then change spaces to _
# then rename alpha/beta/gamma to a/b/g
olink5 <- olink4 %>%
    mutate(Assay = gsub("-", "", Assay)) %>%
    mutate(Assay = gsub(" ", "_", Assay)) %>%
    mutate(Assay = gsub("alpha", "a", Assay)) %>%
    mutate(Assay = gsub("beta", "b", Assay)) %>%
    mutate(Assay = gsub("gamma", "g", Assay))

# count (remove) samples with QC warning
olink5$QC_Warning <- as.factor(olink5$QC_Warning)
str(olink5)
summary(olink5$QC_Warning) #2576 have warning

# olink6 <- olink5 %>%
#     filter(QC_Warning == "Pass")

# change values < LOD of that assay to 50% LOD
olink6 <- olink5 %>%
    mutate(corr_NPX = if_else(LOD > NPX, LOD/2, NPX))


# Create columns for pt_ID, time_point, sample_type for each cohort
# ALCO
ALCO <- olink6 %>%
    filter(cohort == "ALCO") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", sub("-.+$", "", SampleID))) %>%  # substitute 2-7 letters and anything after a "-" with nothing
    mutate(time_point = sub("^.+-T", "", sub("S.+$", "", SampleID))) %>% # add timepoint
    mutate(sample_type = sub("^.+S", "", SampleID)) # add sample_type (peripheral or liver vein)

# Remove duplicate sample
ALCO <- ALCO %>%
    filter(sample_type != "1-d")
# remove samples with QC warning
ALCO <- ALCO %>%
    filter(QC_Warning == "Pass")


# TODO: Add subgroup column: 100=ALD, 200=NAFLD, 300=healthy
# mutate(subgroup = pt_ID < 200))

# Change pt_ID, time_point and sample_type to factor
ALCO$pt_ID <- as.factor(ALCO$pt_ID) # 40 indv
ALCO$time_point <- as.factor(ALCO$time_point) # 3 TPs
ALCO$sample_type <- as.factor(ALCO$sample_type) # 2 sample types


str(ALCO)


# RFX
RFX <- olink5 %>%
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


# TIPS
TIPS <- olink5 %>%
    filter(cohort == "TIPS") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", SampleID)) %>%
    mutate(pt_ID = sub("[A-Z]{2}", "", pt_ID)) %>% #remove the 2 letters
    mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID))













