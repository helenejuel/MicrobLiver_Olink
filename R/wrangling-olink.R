# Load packages
source(here::here("R/package-loading.R"))

# Load data from raw NPX data
olink <- OlinkAnalyze::read_NPX(here::here("data-raw/20202249_Juel_NPX.xlsx"))

# Let the wrangling begin!
# str(olink)

# Fix issue with TIPS naming
# replace typo TIEP80EC	--> TIPS80EC
olink2 <- olink %>%
    mutate(SampleID = sub("TIEP80EC", "TIPS80EC", SampleID))
# replace typo TIPSKC/93 --> TIPS93
olink2 <- olink2 %>%
    mutate(SampleID = sub("TIPSKC/93", "TIPS93", SampleID))

# Add columns for cohort and pt_ID to be able to extract each cohort
# I have a sampleID column with values in the format "cohort-patientID-timepoint-sampletype" that I ultimately need to break down to 4 separate columns.
# The cohort is 2-7 UPPERCASE letters, which I believe is written like this "^[A-Z]{2,7}"
# No separator, then patientID is 1-4 numbers, written like this: "[1-9]{1-4}"
# Only some samples have timepoint and sampletype, and they have different patterns depending on the cohort, so I think the easiest will be to first extract cohort and patientID, and then extract the last 2 for each cohort separately.
#grep selects all rows that have this pattern in the SampleID
#grep("^[A-Z]{2,7}", olink2$SampleID) # from the beginning, 2-7 uppercase letters for cohort
#grep("[1-9]{1,4}", olink2$SampleID) # 1-4 numbers for pt_ID
olink2 <- olink2 %>%
    mutate(cohort = sub("[0-9]{1,4}.*$", "", SampleID)) %>% # substitute 1-4 numbers and everything after this (.*$) with nothing
    mutate(cohort = sub("_.*$", "", cohort)) # substitute anything after a _ with nothing to remove long control names and just keep "CONTROL"

#Keep in mind that * means "match at least 0 times" while + means "match at least 1 time". For some reason the code using .+$ did not replace the single number in the ppts with a single cifre ID.

# Transform to factor to check correct number of cohorts
olink3 <- olink2
olink3$cohort <- as.factor(olink3$cohort)

str(olink3)
summary(olink3$cohort) # if transformed to factor for check: 9 different cohorts, including "Bridge" and "CONTROL"

# Mark duplicate samples in a new column
# olink5 <- olink4 %>%
#     mutate(duplicate = grep("d.+$", SampleID)) #not working
### Do this separately for each cohort

# change cytokine names
# first remove all -
# then change spaces to _
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

# ALCO
# Create columns for pt_ID, time_point, sample_type for each cohort
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


# RFX
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


# TIPS
TIPS <- olink4 %>%
    filter(cohort == "TIPS") %>%
    mutate(pt_ID = sub("^[A-Z]{2,7}", "", SampleID)) %>%
    mutate(pt_ID = sub("[A-Z]{2}", "", pt_ID)) %>% #remove the 2 letters
    mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID))




# GALA --------------------------------------------------------------------
# Load dataset
load(here::here("data/olink4.rda"))

GALA_olink <- olink4 %>%
    filter(cohort == "ALD" | cohort == "HP") %>%
    mutate(sample_type = sub("^.+[0-9]{1,3}", "", SampleID)) # duplicate sample marked as "-d"

# str(as.factor(GALA_olink$sample_type)) # check we have the duplicate samples as "-d"

# Remove duplicate samples
# GALA_olink %>%
#     filter(sample_type == "-d") #276 duplicates (3 samples)
GALA_olink <- GALA_olink %>%
    filter(sample_type != "-d")  # goes from 50692 to 50416 obs

GALA_olink$sample_type <- NULL # remove sample_type column

# remove samples with QC warning
fail <- GALA_olink %>%
    filter(QC_Warning != "Pass") # 1012 obs (12 samples)
unique(fail$SampleID)

GALA_olink <- GALA_olink %>%
    filter(QC_Warning == "Pass") # goes from 50416 to 49404 obs

# Load phenotype data
GALA_SIP_pheno <- read_excel(here::here("data-raw/ALD_HP_SIP_phenotypes.xlsx"), sheet = "small")
str(GALA_SIP_pheno)
# Change gender, abstinent to factor
# for(i in c(2,4)) {
#     GALA_SIP_pheno[,i] = as.factor(GALA_SIP_pheno[,i]} # not working, but should as per old code:
# for(i in c(4:ncol(CHLM))) {CHLM[,i]= as.numeric(CHLM[,i])}
# GALA_SIP_pheno[, "gender"] = as.factor(GALA_SIP_pheno[, "gender"] # also not working
GALA_SIP_pheno$gender <- as.factor(GALA_SIP_pheno$gender) # works
GALA_SIP_pheno$abstinent <- as.factor(GALA_SIP_pheno$abstinent) # also works

#Change age GALA_SIP_pheno[,3] and clinical parameters GALA_SIP_pheno[,c(5:28)] to numeric
GALA_SIP_pheno$age <- as.numeric(GALA_SIP_pheno$age)
# for(j in c(5:28)) {
#     GALA_SIP_pheno[, j] <- as.numeric(GALA_SIP_pheno[, j])
# } # not working, skip for now

# filter for ALD and HP
GALA_pheno <- GALA_SIP_pheno %>%
    mutate(cohort = sub("[0-9]{1,5}", "", ID)) %>%  # substitute 1-5 numbers with nothing
    filter(cohort == "ALD" | cohort == "HP") %>%
    mutate(SampleID = ID) %>% # add SampleID column to be able to merge with Olink data
    mutate(pt_ID = sub("[A-Z]{2,3}", "", ID)) #remove the 2-3 letters in a new pt_ID column to merge with SIP

# str(as.factor(GALA_pheno$cohort)) # check: 2 levels left

# Add SIP samples that included as ALD
SIP <- GALA_SIP_pheno %>%
    mutate(cohort = sub("[0-9]{1,5}", "", ID)) %>%
    filter(cohort == "SIP") %>%
    mutate(pt_ID = sub("[A-Z]{3}", "", ID)) %>% #remove the 3 letters
    mutate(SampleID = paste0("ALD", pt_ID))

# rbind GALA_pheno with SIP
GALA_SIP_pheno2 <- rbind(GALA_pheno, SIP) # back to 1187 individuals
GALA_SIP_pheno2$cohort <- NULL # remove cohort column before merging

# match and merge Olink with phenotype data
GALA <- merge(GALA_olink, GALA_SIP_pheno2, by = "SampleID") # goes from 49404 to 44712 lines from GALA_olink = 51 ppts

GALA_olink_IDs <- sort(unique(GALA_olink$SampleID)) #537 samples
GALA_pheno_IDs <- sort(unique(GALA_pheno$SampleID)) #458 samples --> SIPHON-ALD samples removed from phenotype data! Attempted fixed by matching ID number between GALA_olink and GALA_pheno for SIP only, but still some phenotypes missing or not matched

GALA_IDs <- sort(unique(GALA$SampleID)) #486 samples
TODO: check that HP047 is matched with phenotype data and that there are no other instances of this issue.
TODO: replace phenotype dataset with updated version

str(GALA)

# Graph by Kleiner
GALA %>%
    filter(!is.na(as.numeric(kleinerFscore))) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = kleinerFscore, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("IL8")

# Graph by TE
GALA %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = as.numeric(te), y = corr_NPX)) +
    geom_point(aes(color = cohort)) +
    geom_smooth(method = lm, color = "black") +
    scale_x_log10() +
    ggtitle("IL8")

# Run linear mixed model
GALA$te <- as.numeric(GALA$te)

lm_GALA_te <- olink_lmer(GALA,
           variable = "te",
           outcome = "corr_NPX",
           random = "cohort",
           return.covariates = TRUE)

# Extract list of significant cytokines
sign_cytokines <- lm_GALA_te %>%
    filter(Adjusted_pval < 0.0001) %>%
    pull(Assay)

# Plotting ----------------------------------------------------------------

theme_set(theme_bw())

# histogram
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = corr_NPX)) +
    geom_density() +
    facet_grid(rows = vars(subgroup))

# boxplot
# time points grouped
# Loop over all cytokines
allcytokines_subgroup <- list()

for(i in cytokines[1:92]) {
    ploti <- ALCO %>%
        filter(Assay == i) %>%
        ggplot(aes(x = subgroup, y = corr_NPX)) +
        geom_boxplot() +
        ggtitle(paste(i))

    plotlist[[i]] <- ploti
    # print(plot1)
}

pdf(here::here("doc/images/ALCO_allcytokines.pdf"), height = 60, width = 40) # height = 5 for each row of 8 -> 60 for 12x8
do.call('grid.arrange', c(plotlist, ncol=8))
dev.off()

# plot significant cytokines
# run ANOVA
anova_ALCO_group <- olink_anova(ALCO, variable = "subgroup")

# Extract list of significant cytokines
sign_cytokines <- anova_ALCO_group %>%
    filter(Adjusted_pval < 0.0001) %>%
    pull(Assay)
# print(sign_cytokines) # 51 cytokines significant

signcytokines_subgroup <- list()

for(i in sign_cytokines[1:51]) {
    ploti <- ALCO %>%
        filter(Assay == i) %>%
        ggplot(aes(x = subgroup, y = corr_NPX)) +
        geom_violin() +
        ggtitle(paste(i))

    plotlist2[[i]] <- ploti
    # print(plot1)
}

pdf(here::here("doc/images/ALCO_sign_cytokines_subgroup.pdf"), height = 35, width = 40) # height = 5 for each row of 8 -> 35 for 7x8
do.call('grid.arrange', c(plotlist2, ncol=8))
dev.off()

# Individual graphs
# By subgroup
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX)) +
    geom_boxplot(outlier.colour = "transparent") +
    geom_jitter() +
    ggtitle("IL8")

# time points separated
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX, color = time_point)) +
    geom_boxplot() +
    ggtitle("IL8")

# sample types separated
ALCO %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = subgroup, y = corr_NPX, color = sample_type)) +
    geom_boxplot() +
    ggtitle("IL8")

# only baseline
ALCO %>%
    filter(Assay == "IL8") %>%
    filter(time_point == 1) %>%
    ggplot(aes(x = subgroup, y = corr_NPX)) +
    geom_boxplot(outlier.colour = "transparent") +
    geom_jitter() +
    ggtitle("IL8_baseline")

# Plot using Olink package
# plot all cytokines
#ALCO_subgroup_1 <-
olink_boxplot(ALCO,
              variable = "subgroup",
              olinkid_list = cytokine_IDs[1:16],
              verbose = T,
              number_of_proteins_per_plot = 16)
# ggsave(here::here("doc/images/ALCO_subgroup_1.pdf"),
       # ALCO_subgroup_1, width = 10, height = 6) # not working because it is not a ggplot

# plot significant cytokines
anova_ALCO_group <- olink_anova(ALCO, variable = "subgroup")
sign_cytokines <- anova_ALCO_group %>%
    filter(Adjusted_pval < 0.001) %>%
    pull(OlinkID)

olink_boxplot(ALCO,
              variable = "subgroup",
              olinkid_list = sign_cytokines,
              verbose = T,
              number_of_proteins_per_plot = 6) # only prints the last set of plots

