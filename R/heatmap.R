# Load packages
source(here::here("R/package-loading.R"))

# Heatmap -----------------------------------------------------------------

# Aim: create heatmap of all 92 cytokines' corr.NPX values against each other

# Load data
load(here::here("data/GALA.rda"))

# Need to have each cytokine in a column
# dput(names(GALA_heatmap)) # get list of column names
GALA_heatmap <- pivot_wider(data = GALA,
                            id_cols = c(SampleID,
                                        cohort, gender, # factor
                                        age, # numeric
                                        kleiner, fibrosis, nas_inflam, nas_steatosis, te_fibrosis_6, te_fibrosis_8, inflam, steatosis_binary,  # liver parameters as factors
                                        te, kleiner_numeric, inflam_numeric, steatosis_numeric, elf, # numeric liver parameters
                                        bmi, hr, map, # clinical numerical variables
                                        iga, igg, igm, # immunoglobulins (numeric)
                                        alt, ast, ggt, # liver enzymes (numeric)
                                        crp, # C-reactive protein (numeric)
                                        ldl, hdl, trigly, # serum lipids (numeric)
                                        homair, # insulin resistance (numeric)
                                        proc3, # numeric
                                        abstinent, overuse, excess_drink_followup, # factors
                                        hospInfection, liverrelated_event, allMortality, any_event_tot, # events during followup (factor)
                                        days_to_hospInf2, days_to_LRE2, days_to_mort2, # numeric, NAs carried over and "no" = 2xmax
                                        hospInf_1yr, LRE_1yr, mort_1yr, any_event_1yr, # events during 1st year (factor)
                                        LYPLAL1, MARC1, ERLIN1, GPAM, SERPINA1, APOH, TM6SF2, APOC1, MBOAT7, GCKR, PNPLA3, PPARG, HSD17B13, MTTP, SLC39A8, TRIB1, AKNA, # SNPs (numeric)
                                        PRS_alcohol, PRS_BMI, PRS_CAD, PRS_CHD, PRS_CRP, PRS_GGT, PRS_HDL, PRS_height, PRS_IBD, PRS_IL6, PRS_LDL, PRS_ALT, PRS_risk, PRS_stroke, PRS_T2D, PRS_TC, PRS_TG), # PRSs (numeric)
                            # id_cols chooses which columns to keep, in addition to the cytokines
                            names_from = Assay,
                            values_from = NPX)
# Results in df with 531 obs (number of participants) and 92 variables + the number of variables added as id_cols = 172

# str(GALA_heatmap[ , c(1:10)]) # check

# Save as df
# GALA_wide <- GALA_heatmap
# usethis::use_data(GALA_wide, overwrite = T)


# Save colnames to be included in heatmap as vector
# all_variables <- names(GALA_heatmap[, 3:ncol(GALA_heatmap)])
all_numeric_variables <- c(all_clinical_numeric, all_genetics, all_cytokines) #150 numeric variables
all_cytokines <- names(GALA_heatmap[, 81:ncol(GALA_heatmap)])
all_clinical_numeric <- names(GALA_heatmap[ , c("age", "te", "kleiner_numeric", "inflam_numeric", "steatosis_numeric", "elf", "bmi", "hr", "map", "iga", "igg", "igm", "alt", "ast", "ggt", "crp", "ldl", "hdl", "trigly", "homair", "proc3", "days_to_hospInf2", "days_to_LRE2", "days_to_mort2")])
all_genetics <- names(GALA_heatmap[ , c("LYPLAL1", "MARC1", "ERLIN1", "GPAM", "SERPINA1", "APOH", "TM6SF2", "APOC1", "MBOAT7", "GCKR", "PNPLA3", "PPARG", "HSD17B13", "MTTP", "SLC39A8", "TRIB1", "AKNA", "PRS_alcohol", "PRS_BMI", "PRS_CAD", "PRS_CHD", "PRS_CRP", "PRS_GGT", "PRS_HDL", "PRS_height", "PRS_IBD", "PRS_IL6", "PRS_LDL", "PRS_ALT", "PRS_risk", "PRS_stroke", "PRS_T2D", "PRS_TC", "PRS_TG")])
clinliver_SNPs <- names(GALA_heatmap[ , c("te", "kleiner_numeric", "inflam_numeric", "steatosis_numeric", "LYPLAL1", "MARC1", "ERLIN1", "GPAM", "SERPINA1", "APOH", "TM6SF2", "APOC1", "MBOAT7", "GCKR", "PNPLA3", "PPARG", "HSD17B13", "MTTP", "SLC39A8", "TRIB1", "AKNA")])
clin_PRSs <- names(GALA_heatmap[ , c("te", "kleiner_numeric", "inflam_numeric", "steatosis_numeric", "bmi", "PRS_BMI", "homair", "PRS_T2D", "hdl", "PRS_HDL", "ldl", "PRS_LDL","PRS_TC", "trigly", "PRS_TG", "PRS_CAD", "PRS_CHD", "PRS_stroke", "alt", "PRS_ALT", "ast", "ggt", "PRS_GGT", "iga", "igg", "igm", "crp", "PRS_CRP", "IL6", "PRS_IL6", "PRS_IBD", "PRS_alcohol", "PRS_risk", "PRS_height")])
input_var_factor <- names(GALA_wide[, c("cohort", "gender", "abstinent", "overuse", "excess_drink_followup")])
outcome_var_numeric <- names(GALA_wide[, c("te", "kleiner_numeric", "inflam_numeric", "steatosis_numeric", "elf", "days_to_hospInf2", "days_to_LRE2", "days_to_mort2")])
outcome_var_factor <- names(GALA_wide[, c("kleiner", "fibrosis", "nas_inflam", "nas_steatosis", "te_fibrosis_6", "te_fibrosis_8", "inflam", "steatosis_binary", "hospInfection", "liverrelated_event", "allMortality", "any_event_tot", "hospInf_1yr", "LRE_1yr", "mort_1yr", "any_event_1yr")])


# Assign selected variables to the vector "variables", which will be used for correlation matrix and plotting
variables <- all_numeric_variables
# variables <- c(all_clinical_numeric, all_genetics)
# variables <- all_cytokines
# variables <- clin_PRSs
# variables <- clinliver_SNPs

# Create correlation matrix
cor_matrix <- round(cor(GALA_heatmap[, variables],
                        method = "spearman",
                        use = "pairwise.complete.obs"),
                    2)
# Remove duplicate values if you want a correlation triangle
# cor_matrix[lower.tri(cor_matrix)] <- NA

# Transform to a df with 3 columns: Var1, Var2, value (= correlation R2)
# library(reshape2)
cor_matrix <- melt(cor_matrix)

# Change variable names to characters
cor_matrix$Var1 <- as.character(cor_matrix$Var1)
cor_matrix$Var2 <- as.character(cor_matrix$Var2)

# Remove NAs (duplicates) if needed
# cor_matrix <- na.omit(cor_matrix)
# str(cor_matrix) # check

# Full heatmap (too many to really see)
cor_matrix %>%
    ggplot(aes(Var2, Var1)) +
    geom_tile(aes(fill=value), color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint=0, limit=c(-1,1), name="Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=15, hjust=1),
          axis.text.y = element_text(size=15),
          axis.title = element_blank()) +
    ggtitle("Clinical liver variables + SNPs") +
    ylim(variables2) + xlim(variables2) + # Select variables vector for xlim/ylim, to ensure a nice triangular heatmap
    coord_equal()

# Plot, then save
# ggsave(here::here("doc/images/heatmap_all_numeric.jpg"), height = 12, width = 14)
# ggsave(here::here("doc/images/heatmap_clin_PRS.jpg"), height = 12, width = 14)
ggsave(here::here("doc/images/heatmap_clin_SNP.jpg"), height = 12, width = 14)

# Order cytokines in the order highest to lowest R with HGF
HGF <- cor_matrix %>%
    filter(Var1 == "HGF" | Var2 == "HGF") %>%
    mutate(Var3 = ifelse(Var1 == "HGF", Var2, Var1)) %>%
    arrange(desc(value))
variables2 <- unique(HGF$Var3)
# setdiff(variables, variables2) # check that we have included all cytokines

# Order cytokines in the order highest to lowest R with te
te <- cor_matrix %>%
    filter(Var1 == "te" | Var2 == "te") %>%
    mutate(Var3 = ifelse(Var1 == "te", Var2, Var1)) %>%
    arrange(desc(value))
variables2 <- unique(te$Var3)

# Plot, then save
# ggsave(here::here("doc/images/all_cyto_ordered.jpg"),
#                   height = 12, width = 14)
# ggsave(here::here("doc/images/all_bloodmarkers_ordered.jpg"),
#                   height = 12, width = 14)


### Graph only certain correlations
# Extract all correlations with te with value > 0.5
predictors <- cor_matrix %>% # Requires cor_matrix to include both lower and upper tri
    filter(Var2 == "days_to_LRE2" | Var1 == "days_to_LRE2") %>%
    filter(value > 0.4 | value < -0.4) %>%
    arrange(desc(value)) #use desc(value) for descending
variables3 <- unique(predictors$Var1)

cor_matrix %>%
    ggplot(aes(Var2, Var1)) +
    geom_tile(aes(fill=value), color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint=0, limit=c(-1,1), name="Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=1),
          axis.text.y = element_text(size=10),
          axis.title = element_blank()) +
    ggtitle("Correlation with time to liver-related event R>|0.4|") +
    ylim(variables3) + xlim(variables3) +
    coord_equal()


### Look at cytokines correlated with te or other phenotype
# Create correlation matrix as above, including all numeric variables, DO NOT remove duplicate values (lower tri)
# extract all rows containing phenotype of interest from the correlation matrix
# bmi_cor <- cor_matrix %>%
#     filter(Var1 == "bmi") # only homair correlates >0.3
te_cor <- cor_matrix %>%
    filter(Var1 == "te")
hospInf_cor <- cor_matrix %>%
    filter(Var1 == "days_to_hospInf") %>%
    arrange(value) # nothing correlates with days to hospInf >|0.3|, but steatosis = -0.29 and FGF21 = -0.21
steatosis_cor <- cor_matrix %>%
    filter(Var1 == "nas_steatosis") # 17 variables have R>|0.3|
inflam_cor <- cor_matrix %>%
    filter(Var1 == "nas_inflam") # 25 variables have R>|0.4|

# Extract all values above a set threshold
predictors <- inflam_cor %>%
    filter(value > 0.4  | value < -0.4) %>%
    arrange(desc(value))
predictors # check

# To copy table
# write.table(predictors,"clipboard", sep="\t")

# Assign which variables to include in graph
sign_predictors <- predictors$Var2 # significant variables
all_numeric_variables
all_cytokines
phenotypes <- names(GALA_heatmap[, c(4,5,9,6,8,7,10, 26)]) #for TE plot
phenotypes <- names(GALA_heatmap[, c(7,6,4,5,9,8,10, 26)]) #for steatosis plot
phenotypes <- names(GALA_heatmap[, c(6,4,5,9,7,8,10, 26)]) #for inflam plot

# plot all or significant parameters against one parameter
theme_set(theme_bw(base_size = 10))
cor_matrix %>%
    ggplot(aes(Var2, Var1)) + #change data= to relevant correlation
    geom_tile(aes(fill = value), color="white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), name = "Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=1),
          axis.text.y = element_text(size=10),
          axis.title = element_blank()) +
    xlim(sign_predictors) +
    ylim(phenotypes) +
    ggtitle("Variables with correlation R>|0.4| with liver inflammation")

# ggsave(here::here("doc/images/sign_predictors_te.jpg"),
#                   height = 5, width = 12)
# ggsave(here::here("doc/images/sign_predictors_steatosis.jpg"),
#        height = 5, width = 12)
# ggsave(here::here("doc/images/sign_predictors_inflam.jpg"),
#        height = 5, width = 12)


### compare different populations, e.g. low vs high te, ALD vs HP, overuse, hospInfection

low_fibrosis_heatmap <- GALA_heatmap %>%
    filter(te_fibrosis == "low") # 293 individuals
high_fibrosis_heatmap <- GALA_heatmap %>%
    filter(te_fibrosis == "high") # 225 individuals

HP_heatmap <- GALA_heatmap %>%
    filter(cohort == "HP") # 132 individuals
ALD_heatmap <- GALA_heatmap %>%
    filter(cohort == "ALD") # 404 individuals

overuse_heatmap <- ALD_heatmap %>%
    filter(overuse == "yes") # 96 individuals
no_overuse_heatmap <- ALD_heatmap %>%
    filter(overuse == "no") # 233 individuals
hospInf_heatmap <- ALD_heatmap %>%
    filter(hospInfection == "yes") # 121 individuals
no_hospInf_heatmap <- ALD_heatmap %>%
    filter(hospInfection == "no") # 282 individuals

# Create correlation matrix
variables <- all_numeric_variables
cor_matrix <- round(cor(no_hospInf_heatmap[, variables],
                        method = "spearman",
                        use = "pairwise.complete.obs"),
                    2)
cor_matrix <- melt(cor_matrix)
cor_matrix$Var1 <- as.character(cor_matrix$Var1)
cor_matrix$Var2 <- as.character(cor_matrix$Var2)

# Extract parameters for all values with R>|0.4| for te - for the full cohort
predictors <- te_cor %>%
    filter(value > 0.4  | value < -0.4) %>%
    arrange(desc(value))
# Assign which variables to include in graph
sign_predictors <- predictors$Var2 # significant variables (x axis)
phenotypes <- names(GALA_heatmap[, c(4,5,9,6,8,7,10, 26)]) # y axis for TE plot

# plot
cor_matrix %>%
    ggplot(aes(Var2, Var1)) + #change data= to relevant correlation
    geom_tile(aes(fill = value), color="white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), name = "Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=1),
          axis.text.y = element_text(size=10),
          axis.title = element_blank()) +
    xlim(sign_predictors) +
    ylim(phenotypes) +
    ggtitle("Variables with correlation R>|0.4| with TE, for follow-up no hospInf")

# ggsave(here::here("doc/images/predictors_te_lowfibrosis.jpg"),
#                   height = 4, width = 12)
# ggsave(here::here("doc/images/predictors_te_highfibrosis.jpg"),
#        height = 4, width = 12)
# ggsave(here::here("doc/images/predictors_te_overuse.jpg"),
#        height = 4, width = 12)
# ggsave(here::here("doc/images/predictors_te_hospinf.jpg"),
#        height = 4, width = 12)
# ggsave(here::here("doc/images/predictors_te_nohospinf.jpg"),
#        height = 4, width = 12)

# Correlation test
x1 = GALA_heatmap$te
y1 = GALA_heatmap$IL8
cor.test(x1, y1, alternative = c("two.sided"), method = c("spearman"), exact = NULL)

# Loop plots in base R with log10 axes and lm
# for(i in peptides) {
#     x1 = orderedData2[,"Neut.D126"] #change to functional assay of interest
#     y1 = orderedData2[,i]
#     mod1 <- lm(log10(y1) ~ log10(x1))
#     cor1 <- cor.test(x1, y1, alternative=c("two.sided"), method=c("spearman"), exact=NULL)
#     plot(x=x1, y=y1, log="xy", xlab="", ylab="", main=paste0(i, " rho=", round(cor1$estimate["rho"],2)))
#     abline(mod1)
# }

# add p-values within the main
# "p=",round(cor1$p.value, 4)




