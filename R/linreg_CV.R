# Load packages
source(here::here("R/package-loading.R"))

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
