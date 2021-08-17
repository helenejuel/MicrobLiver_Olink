# Load packages
source(here::here("R/package-loading.R"))

# PCA analysis ------------------------------------------------------------

# Load data
load(here::here("data/GALA_wide.rda"))
str(GALA_wide)

# PCA_exclude <- which(is.na(Norm.test)) # to remove problematic columns, e.g. all values are identical or not numeric
PCA_exclude <- c(1:3, 23:25)# to remove ID, cohort, categorical variables
PCA_exclude <- c(1:26) # to remove all except cytokines

# t(names(GALA_wide)) # to see column name numbers

PCA_data <- GALA_wide[, -PCA_exclude] #All numerical attributes

### Standardization
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
pca.data <- prcomp(standard.corrected.data, center = F, scale. = F) # center and scale already performed manually
summary(pca.data) # to view cumulative proportion of variance explained
pca.data # to view all principal component directions
PC.Directions <- as.data.frame(pca.data$rotation)
PC.Directions[, 1:4] %>%
    arrange(desc(PC1)) #view weights of first 4 PCs
# library(openxlsx)
write.xlsx(pca.data, file="GALA_cytokines_PCAdirections.xlsx", row.names = F)

qplot(x = 1:dim(PCA_data)[2],
      y = summary(pca.data)$importance[3,],
      ylab = "Cumulative variance",
      xlab = "Number of PCs included",
      main = "PCA cytokines")

# Plot first two principal components
qplot(PC1, PC2, data=as.data.frame(pca.data$x))

# Add phenotypes to be able to add labels and colors to the plot
id <- 1:dim(PCA_data)[1]
id <- GALA_wide[, 1]
endpoints <- GALA_wide[, c(3:7, 23:25)] # categorical endpoints te_fibrosis, abstinent, overuse, hospInf, te, kleiner, nas_inflam, nas_steatosis are numerical endpoints
str(endpoints)
plot.pca <- data.frame(pca.data$x, id, endpoints)


# color points according to categorical endpoint
plot.pca %>%
    ggplot(aes(PC1, PC3)) + # try different combinations of PCs
    geom_point(aes(color = te_fibrosis), size = 2) +
    # geom_point(aes(color = hospInfection), size = 2)
    # geom_point(aes(color = overuse), size = 2)
    geom_ellipse(aes(color = te_fibrosis)) # TODO: find correct geom_

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

