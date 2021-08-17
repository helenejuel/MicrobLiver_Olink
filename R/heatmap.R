# Load packages
source(here::here("R/package-loading.R"))

# Heatmap -----------------------------------------------------------------

# Aim: create heatmap of all 92 cytokines' corr.NPX values against each other
# Need to have each cytokine in a column
# Load data
load(here::here("data/GALA.rda"))

GALA_heatmap <- pivot_wider(data = GALA,
                            id_cols = c(SampleID, cohort, te_fibrosis, te, kleiner, nas_inflam, nas_steatosis, meld, elf, bmi, iga, igg, igm, alt, ast, ggt, crp, ldl, hdl, trigly, homair, proc3, abstinent, overuse, hospInfection, days_to_hospInf), # This line chooses which columns to keep, in addition to the cytokines
                            names_from = Assay,
                            values_from = corr_NPX)
# Results in df with 536 obs (number of participants) and 92 variables + the number of variables added as id_cols = 118
# Change kleiner and nas_inflam to numeric
GALA_heatmap <- GALA_heatmap %>%
    mutate(across(c(kleiner, nas_inflam, nas_steatosis), as.numeric))
str(GALA_heatmap[ , c(1:10)]) # check

# Save colnames to be included in heatmap as vector
# all_variables <- names(GALA_heatmap[, 3:ncol(GALA_heatmap)])
all_numeric_variables <- names(GALA_heatmap[, c(4:22, 26, 27:ncol(GALA_heatmap))])
all_blood_markers <- names(GALA_heatmap[ , c(11:22, 27:ncol(GALA_heatmap))])
all_cytokines <- names(GALA_heatmap[, 27:ncol(GALA_heatmap)])
# test_variables <- names(GALA_heatmap[, c(5,10,12,16)]) # select a smaller number of variables

# Assign selected variables to the vector "variables", which will be used for correlation matrix and plotting
variables <- all_numeric_variables
# variables <- all_blood_markers
# variables <- all_cytokines

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
cor_matrix %>%
    ggplot(aes(Var2, Var1)) +
    geom_tile(aes(fill=value), color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint=0, limit=c(-1,1), name="Correlation\n(Spearman)") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=6, hjust=1),
          axis.text.y = element_text(size=6),
          axis.title = element_blank()) +
    ggtitle("All blood markers - ProC3") +
    ylim(variables2) + xlim(variables2) + # Select variables vector for xlim/ylim, to ensure a nice triangular heatmap
    coord_equal()

# Order cytokines in the order highest to lowest R with HGF
HGF <- cor_matrix %>%
    filter(Var1 == "HGF" | Var2 == "HGF") %>%
    mutate(Var3 = ifelse(Var1 == "HGF", Var2, Var1)) %>%
    arrange(desc(value))
variables2 <- unique(HGF$Var3)
# setdiff(variables, variables2) # check that we have included all cytokines

# Order cytokines in the order highest to lowest R with proc3
proc3 <- cor_matrix %>%
    filter(Var1 == "proc3" | Var2 == "proc3") %>%
    mutate(Var3 = ifelse(Var1 == "proc3", Var2, Var1)) %>%
    arrange(desc(value))
variables2 <- unique(proc3$Var3)

# Plot, then save
# ggsave(here::here("doc/images/all_cyto_ordered.jpg"),
#                   height = 12, width = 14)
# ggsave(here::here("doc/images/all_bloodmarkers_ordered.jpg"),
#                   height = 12, width = 14)


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




