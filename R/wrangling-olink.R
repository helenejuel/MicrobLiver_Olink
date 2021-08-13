# Load packages
source(here::here("R/package-loading.R"))


# Load raw data -> olink4 -------------------------------------------------

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


# Alcochallenge -----------------------------------------------------------

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


### Load phenotype data
# TODO: Change to new QC'ed phenotype dataset
GALA_pheno2 <- read_excel(here::here("data-raw/ALD_HP_outcome_HBJ.xlsx"))

# Change data columns to factor
GALA_pheno2 <- GALA_pheno2 %>%
  mutate(across(c(cohort, gender, fibrosis, inflam, abstinent, overuse, alcoholyears, abstinenceyears, liverrelated_event, hospInfection, allMortality, excess_drink_followup), factor)) # consider including kleiner, starts_with("nas")
# TODO: order factors

# Change data columns to numeric
GALA_pheno2 <- GALA_pheno2 %>%
  mutate(across(c(cpa, meld, elf, te, packyears, height, weight, bmi, waist, hip, whr, hr, map, sbp, dbp, alt, ast, alk, bili, chol, crp, ggt, glc, hba1c, hdl, iga, igg, igm, ldl, leu, trigly, insulin, homair, cpeptid, proc3, days_to_LRE, days_to_hospInf, days_to_mort), as.numeric))

# Change column name from CBMR_ID to SampleID to allow merging with olink data
GALA_pheno2 <- GALA_pheno2 %>%
  rename(SampleID = CBMR_ID)

str(GALA_pheno2)

# Match and merge Olink with phenotype data
GALA <- merge(GALA_olink, GALA_pheno2, by = "SampleID") # goes from 49404 to 49312 obs from GALA_olink = 1 ppt!
# Find missing sample
print(setdiff(unique(GALA_olink$SampleID), GALA_pheno2$SampleID)) # ALD1302
# TODO: include phenotypes for ALD1302

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



# Add column for fibrosis level high or low
GALA <- GALA %>%
    mutate(te_fibrosis = if_else(te < 6, "low", "high"))
GALA$te_fibrosis <- as.factor(GALA$te_fibrosis)

# save dataset in data folder
usethis::use_data(GALA, overwrite = T)


# Add outcome data
# Load outcome data
# TODO: Change to new QC'ed phenotype dataset
GALA_outcome <- read_excel(here::here("data-raw/GALA-ALD_outcome_data.xlsx"))
# str(GALA_outcome) # columns types messed up

# Change id to integer
GALA_outcome$id <- as.integer(GALA_outcome$id)
# Change starts_with(days_to) columns to numeric
GALA_outcome <- GALA_outcome %>%
  mutate(across(starts_with("days_to"), as.numeric))
# Change liver-related_event, hospInfection, allMortality, excess_drink_followup columns to factor
GALA_outcome <- GALA_outcome %>%
  mutate(across(c(liverrelated_event, hospInfection, allMortality, excess_drink_followup), factor))
str(GALA_outcome) # check

# Merge GALA with outcome data
# Create SampleID column in outcome data
GALA_outcome <- GALA_outcome %>%
  mutate(SampleID = CBMR_ID)
# Merge
GALA_FU <- merge(GALA, GALA_outcome, by = "SampleID") # 34500 obs = 375 individuals -> dropped from 462 in outcome data and 537 in olink data?
str(GALA_FU)
GALA_FU_IDs <- sort(unique(GALA_FU$SampleID)) #375 samples
# Compare the character vectors
GALA_IDs %in% GALA_FU_IDs
# TODO: extract the FALSE Sample IDs from the comparison
GALA %>%
  filter(SampleID == "ALD2693")

# save dataset in data folder
usethis::use_data(GALA_FU, overwrite = T)


# Graph by Kleiner

# Load dataset
load(here::here("data/GALA.rda"))

GALA %>%
    filter(!is.na(as.numeric(kleinerFscore))) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = kleinerFscore, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
    ggtitle("IL8")

GALA %>%
  filter(Assay == "IL8") %>%
  ggplot(aes(x = kleiner, y = corr_NPX, color = cohort)) +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_jitter(position = position_jitterdodge(), size = 0.5) +
  ggtitle("IL8")

# Graph by TE
GALA %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = as.numeric(te), y = corr_NPX)) +
    geom_point(aes(color = cohort)) +
    geom_smooth(method = lm, color = "black") +
    scale_x_log10() +
    ggtitle("IL8")

# Graph by dichotomized fibrosis
GALA %>%
    filter(!is.na(te_fibrosis)) %>%
    filter(Assay == "IL8") %>%
    ggplot(aes(x = te_fibrosis, y = corr_NPX)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter() +
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


# Heatmap -----------------------------------------------------------------

# Aim: create heatmap of all 92 cytokines' corr.NPX values against each other
# Need to have each cytokine in a column
GALA_heatmap <- pivot_wider(data = GALA,
                            id_cols = c(SampleID, cohort, te, te_fibrosis, kleiner, nas_inflam, inflam, bmi, abstinent, overuse, iga, igg, igm, homair, proc3, hospInfection), # This chooses which columns to keep
                            names_from = Assay,
                            values_from = corr_NPX)
# Results in df with 486 obs (number of participants) and 92 variables + the number of variables added as id_cols = 96
# Change kleiner and nas_inflam to numeric
GALA_heatmap <- GALA_heatmap %>%
  mutate(across(c(kleiner, nas_inflam), as.numeric))
str(GALA_heatmap) # check

# Save colnames to be included in heatmap as vector
# all_variables <- names(GALA_heatmap[, 3:ncol(GALA_heatmap)])
# all_numeric_variables <- names(GALA_heatmap[, c(3,5,6,8, 11:15, 17:ncol(GALA_heatmap))])
all_cytokines <- names(GALA_heatmap[, 17:ncol(GALA_heatmap)])
# test_variables <- names(GALA_heatmap[, c(5,10,12,16)]) # select a smaller number of variables

# Assign selected variables to the vector "variables", which will be used for correlation matrix and plotting
variables <- all_numeric_variables

# Create correlation matrix
cor_matrix <- round(cor(GALA_heatmap[, variables],
                        method = "spearman",
                        use = "pairwise.complete.obs"),
                    2)
# Remove duplicate values if you want a correlation triangle
cor_matrix[lower.tri(cor_matrix)] <- NA

# Transform to a df with 3 columns: Var1, Var2, value (= correlation R2)
# library(reshape2)
cor_matrix <- melt(cor_matrix)
# Change variable names to characters
cor_matrix$Var1 <- as.character(cor_matrix$Var1)
cor_matrix$Var2 <- as.character(cor_matrix$Var2)
# Remove NAs (duplicates) if needed
cor_matrix <- na.omit(cor_matrix)
# str(cor_matrix) # check

# Full heatmap (too many to really see)
ggplot(cor_matrix, aes(Var2, Var1)) +
  geom_tile(aes(fill=value), color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1,1), name="Correlation\n(Spearman)") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4, hjust=1),
        axis.text.y = element_text(size=4),
        axis.title = element_blank()) +
  ggtitle("All cytokines") +
  ylim(variables) + xlim(variables) + # Select variables vector for xlim/ylim, to ensure a nice triangular heatmap
  coord_equal()

### Graph only certain correlations
# Extract all correlations with te with value > 0.5
predictors <- cor_matrix %>% # Requires cor_matrix to include both lower and upper tri
  filter(Var2 == "te" | Var1 == "te") %>%
  filter(value > 0.5 | value < -0.5) %>%
  arrange(desc(value))
variables2 <- unique(predictors$Var1)

cor_matrix %>%
  ggplot(aes(Var2, Var1)) +
  geom_tile(aes(fill=value), color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1,1), name="Correlation\n(Spearman)") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_blank()) +
  ggtitle("Strongest correlations with TE") +
  ylim(variables2) + xlim(variables2) +
  coord_equal()

# TODO: use stats result to assign cytokines to be plotted + order by r
# TODO change plot to piping from cor_matrix, should allow to filter for Var1 == te to only include these correlations.
# TODO: then  create heatmap for different populations, e.g. low vs high te, ALD vs HP, overuse, hospInfection
# sign_variables <-

# Extract all correlations with te with value > 0.5
predictors <- cor_matrix %>% # Requires cor_matrix to include both lower and upper tri
  filter(value > 0.5 | value < -0.5) %>%
  arrange(desc(value))
variables3 <- unique(predictors$Var1)

###at cytokines correlated with te or bmi
# Create correlation matrix
cor_matrix <- round(cor(GALA_heatmap[, variables],
                        method = "spearman",
                        use = "pairwise.complete.obs"),
                    2)
# DO NOT remove duplicate values

# Transform to a df with 3 columns: Var1, Var2, value (= correlation R2)
cor_matrix <- melt(cor_matrix)
# Change variable names to characters
cor_matrix$Var1 <- as.character(cor_matrix$Var1)
cor_matrix$Var2 <- as.character(cor_matrix$Var2)

# extract all rows containing te or bmi from the correlation matrix
bmi_cor <- cor_matrix %>%
    filter(Var1 == "bmi")
te_cor <- cor_matrix %>%
    filter(Var1 == "te")

# Extract all values above a set threshold, removing self-correlation
predictors <- te_cor %>%
    filter(value > 0.5  & value < 1 | value < -0.5) %>%
    arrange(desc(value))
# predictors # check

# To copy table
# write.table(predictors,"clipboard",sep="\t")

# Assign which variables to include in graph
variables <- predictors$Var2 # significant variables
variables <- all_variables

# plot all or significant cytokines against one assay
theme_set(theme_bw(base_size = 8))
ggplot(te_cor, aes(Var2, Var1)) + #change data= to relevant correlation
    geom_tile(aes(fill = value), color="white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), name = "Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4, hjust=1),
          axis.text.y = element_text(size=4), axis.title = element_blank())
    # # ggtitle("All vaccinees") +
    # xlim(variables)
# save as image 1500x100

#graph individual correlations ggplot
# theme_set(theme_bw(base_size = 14))
# ggplot(data=orderedData2, aes(x=Neut, y=orderedData2[,"033VD1"])) +
# geom_point(size=4) +
# geom_smooth(method=lm, color="black") +
# scale_y_log10()

# loop plots
graphics.off()
par(mfrow=c(3,4), mar=c(2,3,2,1), oma=c(0,0,2,0))
# CAFpeptides <- c("033VD1", "034VD1", "036VD1", "180", "185", "191", "261VD4", "262VD4", "263VD4", "264VD4")
peptides <- predictors[c(1:nrow(predictors)-1),"Var1"] #select row 1 to nrow-1 to remove "Phago"
# peptides <- predictors[c(1:nrow(predictors)),"Var1"] #select nrow if self-correlation has already been removed
# simple plots
# for(i in peptides) {
#   plot(x=orderedData4[,"Neut"], y=orderedData4[,i], ylim = c(0,66000), main=paste(peptides[i], "Neut"), xlab="", ylab="")
# }

# log10 axes and lm
for(i in peptides) {
    x1 = orderedData2[,"Neut.D126"] #change to functional assay of interest
    y1 = orderedData2[,i]
    mod1 <- lm(log10(y1) ~ log10(x1))
    cor1 <- cor.test(x1, y1, alternative=c("two.sided"), method=c("spearman"), exact=NULL)
    plot(x=x1, y=y1, log="xy", xlab="", ylab="", main=paste0(i, " rho=", round(cor1$estimate["rho"],2)))
    abline(mod1)
}

# add p-values within the main
# "p=",round(cor1$p.value, 4)
# 1 row: 600x200
# 3 rows: 600x

# Correlation test
x1 = orderedData2[,"Phago"]
y1 = orderedData2[,40] #40=peptide99
cor.test(x1, y1, alternative=c("two.sided"), method=c("spearman"), exact=NULL)

# get sequence of selected peptides - subtract 61 from peptide_number to choose line
annotation[264,c(1:3)]
annotation[c(261,264),c(1:3)] #all vaccinees Neut
annotation[c(254,257,259,261),c(1:3)] #Alum Neut
annotation[c(33,34,36,180,185,191,261:264),c(1:3)] #CAF Neut

write.table(annotation[c(33,34,36,180,185,191,261:264),c(1:3)],"clipboard",sep="\t")


# PCA analysis ------------------------------------------------------------
### NOT WORKING!
# PCA_exclude<-which(is.na(Norm.test)) # to remove problematic columns, e.g. all values are identical or not numeric
PCA_exclude <- c(1:3, 28:99, 121:144)# to remove ID, PS, vaccine, mucosa data except muc.norm2

# t(names(VG))

PCA.data <- VG[, -PCA_exclude] #All numeric attributes, all vaccinated

# Standardization
means <- colMeans(PCA.data, na.rm=T) # find means
datzeromean <- t(apply(PCA.data,1,'-',means)) #subtract means from all values
colMeans(datzeromean,na.rm=T) #colmeans should now be ~0
standard.deviation <- apply(PCA.data,2,sd,na.rm=T) #find SDs
standard.data <- t(apply(datzeromean,1,"/",standard.deviation)) #divide normalised data with SD
standard.corrected.data <- standard.data
apply(standard.corrected.data,2,sd, na.rm = T) #SDs should now be 1

for(i in 1:dim(standard.corrected.data)[2]){
    standard.corrected.data[,i][is.na(standard.data[,i])]<-means[i]/standard.deviation[i] # substitute all NA with means/sd for the specific column
}
sum(is.na(standard.corrected.data)) # Test for any NAs = 0

All.identical.test<-rep(NA,dim(standard.corrected.data)[2]) # Test for any columns with only identical values (if sd = 0 then scale will divide the data with 0 = inf
for(i in 1:dim(standard.corrected.data)[2]){
    All.identical.test[i]<-length(unique(standard.corrected.data[,i]))
}
standard.corrected.data[,which(All.identical.test==1)]
All.identical.names<-colnames(standard.corrected.data)[which(All.identical.test==1)]
# remove columns that cannot be used for pca
if(length(which(All.identical.test==1))==0){
    print("No columns deleted")
} else {
    standard.corrected.data<-standard.corrected.data[,-which(All.identical.test==1)]
    cat(All.identical.names,"has been removed from PCA analysis")
}

# Use prcomp of the standardized corrected data
pca.data<-prcomp(standard.corrected.data, center=F,scale.=F) # center and scale already performed manually
summary(pca.data) # to view cumulative proportion of variance explained
pca.data # to view all principal component directions, copy paste to notepad and saved as txt file
# PC.Directions<-pca.data$rotation
# PC.Directions[,1:4] #view weights of first 4 PCs
# library(openxlsx)
write.xlsx(pca.data, file="2019-03-22 PCA directions.xlsx", row.names = F)

qplot(x=1:dim(PCA.data)[2], y=summary(pca.data)$importance[3,], ylab="Cumulative variance", xlab="Number of PCs included", main="PCA vaccine groups")

# Plot first two principal components
qplot(PC1,PC2,data=as.data.frame(pca.data$x))

# Add labels and colors
kid <- 1:dim(PCA.data)[1]
endpoints <- gamk[,c(2:7)] #all kids
endpoints <- SUSC[,c(2:7)] #only susceptible
str(endpoints)
plot.pca<-data.frame(pca.data$x, kid, endpoints)
# plot.pca<-data.frame(pca.data$x,Animal,endpoints,Clear=ggdata$Clear,IFU.D3=ggdata$IFU.D3)

# Label kid number to points
qplot(PC1,PC3, data=plot.pca, geom=c("point","text"), vjust=(-1), label=kid)

# color points according to endpoint
qplot(PC1,PC2, color=susceptible, data=plot.pca, main="PCA all kids")
ggplot(data=plot.pca, aes(x=PC1, y=PC2)) + geom_point(aes(color=Village), size=4)
qplot(PC1,PC2, color=inf, data=plot.pca)
qplot(PC1,PC2, color=dis, data=plot.pca)
ggplot(data=plot.pca, aes(x=PC1, y=PC2)) +
    geom_point(aes(color=log10(load)), size=4) +
    scale_color_gradient(low="blue", high="red")
ggplot(data=plot.pca, aes(x=PC2, y=PC5)) +
    geom_point(aes(color=inf.lgth.med, shape=Village), size=4) +
    scale_color_gradient(low="blue", high="red")
ggplot(data=plot.pca, aes(x=PC1, y=PC2)) +
    geom_point(aes(fill=inf.lgth.med, size=susceptible, color=Village), shape=21) + #21 has line and fill
    scale_fill_gradient(low="blue", high="red")


# combined plots
p1 <- qplot(PC1,PC2, color=susceptible, data=plot.pca) + theme(legend.position = "none")
p2 <- qplot(PC1,PC3, color=susceptible, data=plot.pca) + theme(legend.position = "none")
p3 <- qplot(PC1,PC4, color=susceptible, data=plot.pca) + theme(legend.position = "none")
p4 <- qplot(PC2,PC3, color=susceptible, data=plot.pca) + theme(legend.position = "none")
p5 <- qplot(PC2,PC4, color=susceptible, data=plot.pca) + theme(legend.position = "none")
p6 <- qplot(PC3,PC4, color=susceptible, data=plot.pca) + theme(legend.position = "none")
# library(cowplot)
plot_grid(p1,p2,p3,p4,p5,p6, nrow=2)

# Plot one PC with one end point
qplot(PC1,log10(load), color=susceptible, data=plot.pca)
qplot(PC1, log10(inf.lgth.med), color=Village, shape=susceptible, data=plot.pca)
qplot(PC1, susceptible, color=Village, data=plot.pca)


# Linear prediction models --------------------------------------------------------------

# select all SvD peptides from the CAF group
orderedData1 <- JPTdata[t(pdata[,"Vaccine"]) %in% c("CAF")]
# orderedData1 <- JPTdata[t(pdata[,"Vaccine"]) %in% c("Alum")]
# orderedData1 <- JPTdata[t(pdata[,"Vaccine"]) %in% c("CAF", "Alum")]
# transpose and add neutralisation titres
orderedData2 <- as.data.frame(cbind(t(orderedData1), pdata[pdata[,"Vaccine"] %in% c("CAF"), "Mean.titre.D126.F"]))
p1 <- c(paste0(annotation[,"CTH522.peptide.no"], annotation[,"Domain"]), "Neut")
names(orderedData2) <- p1
# Find NAs
sum(is.na(CAF)) #2 NAs - both in Neut titre:
for(i in names(CAF)) {
    print(which(is.na(CAF[,i]), arr.ind = T))
}
which(is.na(CAF[,6]),arr.ind=T) #none in D126 titre
sum(is.na(orderedData2)) #no NAs

# Linear regression automatically removes the whole row if there is an NA, so columns with many NAs should be removed, e.g.
# lindata <- orderedData2[,-20]
lindata <- orderedData2
str(lindata)

### Load toolbox
# setwd('C:/Users/HBJU/Documents/R/R_scripts/02450Toolbox_R') #Should change this to source directly
source("Tools/source_tools.R")
path = "Tools/"
sourceDir(path, exceptions=c("source_tools.R"))
library(cvTools)

### Cross-Validation - Linear regression
# This function takes as input a training and a test set.
#  1. It fits a linear model on the training set using lm.
#  2. It estimates the output of the test set using predict.
#  3. It computes the sum of squared errors.
funLinreg <- function(X_train, y_train, X_test, y_test){
    Xr <- data.frame(X_train)
    Xtest <- data.frame(X_test)
    if(dim(as.matrix(X_train))[2]!=0){
        xnam <- paste("X", 1:dim(as.matrix(X_train))[2], sep="")
        colnames(Xr) <- xnam
        colnames(Xtest) <- xnam
        (fmla <- as.formula(paste("y_train ~ ", paste(xnam, collapse= "+"))))
    }else{
        xnam <- 1
        (fmla <- as.formula(paste("y_train ~ ", paste(xnam, collapse= "+"))))
    }
    mod = lm(fmla, data=Xr)
    preds <- predict(mod, newdata = Xtest)
    sum((y_test-preds)^2)
}

### Two-layer cross validation
# Create crossvalidation partition for evaluation
EOI<-"Neut" # End point of interest
X <-lindata[,c(1:ncol(lindata)-1)]
y <-lindata[,EOI]
attributeNames <- names(lindata[,c(1:ncol(lindata)-1)])
N = dim(lindata)[1]
M = dim(X)[2]
K = 14

# set.seed(42) # for reproducibility, not used for leave-one-out
CV <- cvFolds(N, K=K)
# set up vectors that will store sizes of training and test sizes
CV$TrainSize <- c()
CV$TestSize <- c()

# Initialize variables
Features <- matrix(rep(NA, times=K*M), nrow=K)
Error_train <- matrix(rep(NA, times=K), nrow=K)
Error_test <- matrix(rep(NA, times=K), nrow=K)
Error_train_fs <- matrix(rep(NA, times=K), nrow=K)
Error_test_fs <- matrix(rep(NA, times=K), nrow=K)
Error_train_gs <- matrix(rep(NA, times=K), nrow=K)
Error_test_gs <- matrix(rep(NA, times=K), nrow=K)

Cost.threshold <- 0.01

# Start the clock!
ptm <- proc.time()

# For each crossvalidation fold
for(k in 1:K){
    cat('Crossvalidation fold ', k, '/',K,"\n")

    # Extract the training and test set
    X_train <- X[CV$subsets[CV$which!=k], ];
    y_train <- y[CV$subsets[CV$which!=k]];
    X_test <- X[CV$subsets[CV$which==k], ];
    y_test <- y[CV$subsets[CV$which==k]];
    CV$TrainSize[k] <- length(y_train)
    CV$TestSize[k] <- length(y_test)

    # Use 13-fold crossvalidation for sequential feature selection (only 13 since 1 is used for outer loop)
    fsres <- forwardSelection(funLinreg, X_train, y_train, nfeats=NULL, minCostImprovement=Cost.threshold, stoppingCrit = "minCostImprovement", cvK=K-1);

    # Extract selected features from the forward selection routing
    selected.features <- fsres$featsIncluded

    # Save the selected features
    Features[k,] = fsres$binaryFeatsIncluded

    # Fix to select the right number (lowest squared error) of attributes
    Correction<-which(fsres$costs==min(fsres$costs))
    selected.features <- fsres$featsIncluded[1:Correction]
    binary.selected.features<-names(X_train)%in%names(X_train)[selected.features]
    Features[k,]<-binary.selected.features
    # The function forwardSelection is supposed to terminate if squared error increases when using minCostImprovement. It correctly terminates after it see an increase, but still chooses to use the feature that caused the increase. Fixed 28/10-17
    # Correction<-length(fsres$costs) # run without fix

    # Compute squared error without feature subset selection
    Error_train[k] = funLinreg(X_train, y_train, X_train, y_train);
    Error_test[k] = funLinreg(X_train, y_train, X_test, y_test);

    #Compute squared error of guessing as average y
    Error_train_gs[k] = sum((y_train-mean(y_train))^2)
    Error_test_gs[k] = sum((y_test-mean(y_train))^2)

    # Compute squared error with feature subset selection
    Error_train_fs[k] = funLinreg(X_train[,selected.features], y_train, X_train[,selected.features], y_train);
    Error_test_fs[k] = funLinreg(X_train[,selected.features], y_train, X_test[,selected.features], y_test);

    # Show variable selection history
    # mfig(sprintf('(%d) Feature selection',k));
    I = length(fsres$costs[1:Correction]) # Number of iterations
    par(mfrow=c(1,2))
    # Plot error criterion
    plot(fsres$costs[1:Correction], xlab='Iteration', ylab='Squared error (crossvalidation)', main='Value of error criterion');
    # Plot feature selection sequence
    bmplot(attributeNames, 1:I, fsres$binaryFeatsIncludedMatrix[1:Correction,]);
}
# Stop the clock
proc.time() - ptm

# Display results
print('Linear regression without feature selection:')
print(paste('- Training error:', sum(Error_train)/sum(CV$TrainSize)));
print(paste('- Test error:', sum(Error_test)/sum(CV$TestSize)));

print('Linear regression with sequential feature selection:');
print(paste('- Training error:', sum(Error_train_fs)/sum(CV$TrainSize)));
print(paste('- Test error:', sum(Error_test_fs)/sum(CV$TestSize)));

print('Guessing the average of the output:');
print(paste('- Training error:', sum(Error_train_gs)/sum(CV$TrainSize)));
print(paste('- Test error:', sum(Error_test_gs)/sum(CV$TestSize)));

# Show the selected features
# dev.off()
par(mfrow=c(1,1))
bmplot(attributeNames, 1:K, Features, xlab='Crossvalidation fold', ylab='', main='Attributes selected');
# For array data, bmplot is not appropriate (too many attributes)
# Want to extract the selected features from each fold and save in a df
Features #logical matrix with which peptides were selected in each fold
which(Features[1,]==T) #which features in fold 1 were selected

# Prepare matrix
feat <- matrix(nrow=K, ncol=12)
feat[i,] <- which(Features[i,]==T) #not working

### EXTRACT FEATURES - NEEDS WORK
for(i in c(1:K)){
    print(which(Features[i,]==T))
} #kinda works, but prints ascending

# Calculating Test error using best model only and not the average
best.lin.attributes <- c("018", "035VD1", "176", "194", "226")

Error_lin_train<-rep(NA,K)
Error_lin_test<-rep(NA,K)
for(k in 1:K){
    cat('Crossvalidation fold ', k, '/',K,"\n")

    # Extract the training and test set
    X_train <- X[CV$which!=k, best.lin.attributes];
    y_train <- y[CV$which!=k];
    X_test <- X[CV$which==k, best.lin.attributes];
    y_test <- y[CV$which==k];
    CV$TrainSize[k] <- length(y_train)
    CV$TestSize[k] <- length(y_test)

    # Compute squared error
    Error_lin_train[k] = funLinreg(X_train, y_train, X_train, y_train);
    Error_lin_test[k] = funLinreg(X_train, y_train, X_test, y_test);
}

print('Linear regression with 5 best attibutes:');
print(paste('- Training error:', sum(Error_lin_train)/sum(CV$TrainSize)));
print(paste('- Test error:', sum(Error_lin_test)/sum(CV$TestSize)));

Best.lin.test.error<-sum(Error_lin_test)/sum(CV$TestSize)

### Check a model ###
#5-parameter
# best.model <- lm(formula = Neut ~ 018+035VD1+176+194+226, data=lindata) #not working: column names that are numbers are not syntactically valid
best.model <- lm(formula = Neut ~ `018` + `035VD1` + `176` + `194` + `226`, data=lindata) #wrap invalid column names in backticks

# Graph model y= a + bx1 + cx2 + dx3 + ex4 + fx5
# model5 <- best.model$coefficients[1] + best.model$coefficients[2]*lindata[,names(best.model$coefficients)[2]] + best.model$coefficients[3]*lindata[,names(best.model$coefficients)[3]] + best.model$coefficients[4]*lindata[,names(best.model$coefficients)[4]] + best.model$coefficients[5]*lindata[,names(best.model$coefficients)[5]] + best.model$coefficients[6]*lindata[,names(best.model$coefficients)[6]] #not working because of the backticks

# best.lin.attributes <- c("018", "035VD1", "176", "194", "226")
model5 <- best.model$coefficients[1] + best.model$coefficients[2]*lindata[,best.lin.attributes[1]] + best.model$coefficients[3]*lindata[,best.lin.attributes[2]] + best.model$coefficients[4]*lindata[,best.lin.attributes[3]] + best.model$coefficients[5]*lindata[,best.lin.attributes[4]] + best.model$coefficients[6]*lindata[,best.lin.attributes[5]]

# 3-parameter
BLA <- c("018", "035VD1", "194")
best.model <- lm(formula = Neut ~ `018` + `035VD1` + `194`, data=lindata)
model3 <- best.model$coefficients[1] + best.model$coefficients[2]*lindata[,BLA[1]] + best.model$coefficients[3]*lindata[,BLA[2]] + best.model$coefficients[4]*lindata[,BLA[3]]

# 2-parameter
BLA <- c("035VD1", "194")
best.model <- lm(formula = Neut ~ `035VD1` + `194`, data=lindata)
model2 <- best.model$coefficients[1] + best.model$coefficients[2]*lindata[,BLA[1]] + best.model$coefficients[3]*lindata[,BLA[2]]

# 1-parameter
BLA <- "035VD1"
best.model <- lm(formula = Neut ~ `035VD1`, data=lindata)
model1 <- best.model$coefficients[1] + best.model$coefficients[2]*lindata[,BLA[1]]

# theme_set(theme_bw(base_size = 14))
ggplot(data= lindata, aes(x=model2, y=Neut)) +
    geom_point(size=3) +
    geom_smooth(method=lm) +
    xlab("Neut_prediction") +
    #ggtitle("5-parameter model: 018, 035VD1, 176, 194, 226")
    #ggtitle("3-parameter model: 018, 035VD1, 194")
    ggtitle("2-parameter model: 035VD1, 194")
#ggtitle("1-parameter model: 035VD1")
#size 600x400

# graph residual error vs attributes
ggplot(data=lindata, aes(x=Ocu.G.W22, y=(Clear-model2))) + geom_point(size=3) + ylab("Residual error")
ggplot(data=lindata, aes(x=Ser.Max, y=(Clear-model2))) + geom_point(size=3) + ylab("Residual error")
