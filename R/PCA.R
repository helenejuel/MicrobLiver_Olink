# Load packages
source(here::here("R/package-loading.R"))

# PCA analysis ------------------------------------------------------------

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
