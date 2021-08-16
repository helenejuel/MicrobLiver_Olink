# Load packages
source(here::here("R/package-loading.R"))

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
