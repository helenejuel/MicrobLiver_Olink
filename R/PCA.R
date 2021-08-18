# Load packages
source(here::here("R/package-loading.R"))

# PCA analysis ------------------------------------------------------------

# Load data
load(here::here("data/GALA_wide.rda"))
str(GALA_wide)

# PCA_exclude <- which(is.na(Norm.test)) # to remove problematic columns, e.g. all values are identical or not numeric
# PCA_exclude <- c(1:3, 23:25)# to remove ID, cohort, categorical variables
PCA_exclude <- c(1:26) # to remove all except cytokines

# t(names(GALA_wide)) # to see column name numbers

PCA_data <- GALA_wide[, -PCA_exclude] #All numerical attributes
str(PCA_data) # check

# Save the untransformed PCA data, to be able to use the PCA_data as object name in the script without having to change it all the way
PCA_data_UT <- PCA_data

# Need to check that all columns have normal distribution, as PCA assumes normally distributed data in each column.
# Check: http://www.sthda.com/english/wiki/normality-test-in-r
# One method is to visualise data with histograms or qq-plots
PCA_data %>%
    ggplot(aes(IL10RB)) +
    geom_density()
# hist(PCA_data$IL10RB) # for baseR plot
PCA_data %>%
    ggplot(aes(sample = IL10RB)) +
    stat_qq() + stat_qq_line()
# qqline(PCA_data$IL10RB) # for baseR plot

# Another method is to use a normality test, like Shapiro-Wilk's. Note that small n will often pass this test, in which case it is important to perform a visual inspection as well.
shapiro.test(PCA_data$IL10RB) # If p>0.05, data are not significantly different from a normal distribution

# Plot qq plots of all cytokines
qq_data <- pivot_longer(PCA_data, cols = c(1:ncol(PCA_data)), names_to = "Assay")
# ggplot(qq_data, aes(x = value)) +
#     geom_density() +
#     facet_wrap(vars(Assay), nrow = 8) #looks a little weird?
ggplot(qq_data, aes(sample = value)) +
    stat_qq() + stat_qq_line() +
    facet_wrap(vars(Assay), nrow = 8)
ggsave(here::here("doc/images/all_cyto_qq.pdf"), height = 12, width = 14)



# Apply Shapiro test to all 92 cytokines
# save the Shapiro test in a list of lists, with a list for each cytokine
shapiro <- apply(PCA_data, 2, shapiro.test)
# shapiro <- apply(PCA_data_log, 2, shapiro.test)
shapiro <- apply(PCA_data_INT, 2, shapiro.test)
# extract p values from this list of lists using do.call
shapiro_p <- as.data.frame(do.call(rbind, lapply(shapiro, function(v){v$p.value})))
# add column with cytokine names
shapiro_p <- shapiro_p %>%
    mutate(Assay = cytokines)
# Extract cytokines with p.value <0.5
non_normal <- shapiro_p %>%
    filter(V1 < 0.05)
non_normal # view



# Since 87 of the 92 cytokines were not normally distributed, I will apply log10 transformation and run the shapiro test again
PCA_data_log <- log10(PCA_data[, c(1:ncol(PCA_data))])

# Still 81 cytokines are not normally distributed
# I will try to apply inverse normal transformation instead
# Test on a single column
# IL10B_INT <- qnorm((rank(PCA_data$IL10RB, na.last = "keep") - 0.5) / sum(!is.na(PCA_data$IL10RB)))
# hist(IL10B_INT) # seems to work
# hist(PCA_data$IL13)
# hist(qnorm((rank(PCA_data$IL13, na.last = "keep") - 0.5) / sum(!is.na(PCA_data$IL13)))) # however, does not seem to perform too well with the cytokines with many measurements below LoD

PCA_data_INT <- PCA_data
for(i in c(1:ncol(PCA_data))) {
    PCA_data_INT[, i] <- qnorm((rank(PCA_data[, i], na.last = "keep") - 0.5) / sum(!is.na(PCA_data[, i])))
} # Also seems to work
# Now only 33 are not normally distributed - seems to be those with many values below LoD

# Will remove these 33 cytokines for now
# Run shapiro test and filter away those cytokines that end in the non_normal$Assay column
PCA_data_INT_clean <- PCA_data_INT[, setdiff(names(PCA_data_INT), non_normal$Assay)]
summary(PCA_data_INT_clean) # 59 cytokines remain - means are already 0
sd(PCA_data_INT_clean$IL10RB) # and sd is ~1
# So no need to standardize the data before PCA


### Standardization
# Choose the dataset to be used in the PCA
PCA_data <- PCA_data_UT

means <- colMeans(PCA_data, na.rm = T) # find means
datzeromean <- t(apply(PCA_data, 1, '-', means)) #subtract means from all values
# colMeans(datzeromean, na.rm = T) # Check: colmeans should now be ~0
standard.deviation <- apply(PCA_data, 2, sd, na.rm = T) #find SDs
standard.data <- t(apply(datzeromean, 1, "/", standard.deviation)) #divide normalised data with SD
standard.corrected.data <- standard.data
# apply(standard.corrected.data, 2, sd, na.rm = T) #Check: SDs should now be 1

for(i in 1:dim(standard.corrected.data)[2]) {
    standard.corrected.data[,i][is.na(standard.data[,i])] <- means[i]/standard.deviation[i] # substitute all NA with mean/sd for the specific column
    }
# sum(is.na(standard.corrected.data)) # Check: Test that all NAs = 0

# Test for any columns with only identical values
# if sd = 0 then scale will divide the data with 0 = inf = error
# Prepare list
All.identical.test <- rep(NA, dim(standard.corrected.data)[2])
# Find columns that have only 1 unique datapoint
for(i in 1:dim(standard.corrected.data)[2]){
    All.identical.test[i] <- length(unique(standard.corrected.data[,i]))
    }
# list of colnames of assays with identical data
All.identical.names <- colnames(standard.corrected.data)[which(All.identical.test == 1)]
# remove columns that cannot be used for pca
if(length(which(All.identical.test == 1)) == 0) {
    print("No columns deleted")
    } else {
    standard.corrected.data <- standard.corrected.data[, -which(All.identical.test == 1)]
    cat(All.identical.names,"has been removed from PCA analysis")
    }

# Use prcomp on the standardized corrected data
# Assign the INT_clean data directly
standard.corrected.data <- PCA_data_INT_clean
pca.data <- prcomp(standard.corrected.data, center = F, scale. = F) # center and scale already performed manually
summary(pca.data) # to view cumulative proportion of variance explained
pca.data # to view all principal component directions
PC.Directions <- as.data.frame(pca.data$rotation)
PC.Directions[, 1:4] %>%
    arrange(desc(PC1)) #view weights of first 4 PCs
# library(openxlsx)
# write.xlsx(pca.data, file="GALA_cytokines_PCAdirections.xlsx", row.names = F)

qplot(x = 1:dim(standard.corrected.data)[2],
      y = summary(pca.data)$importance[3,],
      ylab = "Cumulative variance",
      xlab = "Number of PCs included",
      main = "PCA cytokines")

# Plot first two principal components
qplot(PC1, PC2, data=as.data.frame(pca.data$x))

# Add phenotypes to be able to add labels and colors to the plot
id <- GALA_wide[, 1]
endpoints <- GALA_wide[, c(3:7, 23:25)] # categorical endpoints te_fibrosis, abstinent, overuse, hospInf, te, kleiner, nas_inflam, nas_steatosis are numerical endpoints
# str(endpoints) # check
plot.pca <- data.frame(pca.data$x, id, endpoints)


# color points according to categorical endpoint
plot.pca$te_fibrosis
plot.pca %>%
    na.omit(te_fibrosis) %>%
    ggplot(aes(PC1, PC3, color = te_fibrosis)) + # try different combinations of PCs
    geom_point(size = 2) +
    # geom_point(aes(color = hospInfection), size = 2)
    # geom_point(aes(color = overuse), size = 2)
    geom_density2d()
# TODO: replace geom_density with simple ellipses

# color points according to continuous endpoint
plot.pca %>%
    ggplot(aes(PC1, PC3)) +
    geom_point(aes(color = te), size = 2) +
    scale_color_gradient(low = "blue", high = "red")

# facet_wrap on categorical endpoint
plot.pca %>%
    # na.omit(abstinent) %>% # omitting all with NA in abstinent, except the 2 with high fibrosis?!?
    ggplot(aes(PC1, PC3)) +
    geom_point(aes(color = te_fibrosis), size = 2) +
    # facet_wrap("abstinent") +
    # facet_wrap("overuse") +
    facet_wrap("hospInfection") +
    ggtitle("Infection requiring hospitalisation during follow-up")
    # ggtitle("Overuse")
    # ggtitle("Abstinent")


# Show 2 endpoints
ggplot(plot.pca, aes(PC1, PC3)) +
    geom_point(aes(color = te, shape = overuse), size=2) +
    scale_color_gradient(low="blue", high="red")

# Show 3 endpoints - not super easy to see
ggplot(plot.pca, aes(PC1, PC2)) +
    geom_point(aes(fill = te, size = abstinent, color = overuse), shape=21) + #21 has line and fill
    scale_fill_gradient(low="blue", high="red")


# All combinations of PC1-4
# p1 <- qplot(PC1,PC2, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# p2 <- qplot(PC1,PC3, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# p3 <- qplot(PC1,PC4, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# p4 <- qplot(PC2,PC3, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# p5 <- qplot(PC2,PC4, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# p6 <- qplot(PC3,PC4, color=te_fibrosis, data=plot.pca) + theme(legend.position = "none")
# # library(cowplot)
# plot_grid(p1,p2,p3,p4,p5,p6, nrow=2)

# Plot one PC with one end point
plot.pca %>%
    ggplot(aes(log10(te), PC1)) +
    geom_point()

