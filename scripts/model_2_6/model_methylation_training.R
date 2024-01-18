#---
#title: "Model 2 training - omic methylation"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

dim(omic_methylation_training)

## remove features with no variance or all missing values

feature.var <- apply(omic_methylation_training,1,var,na.rm=T)
omic_methylation_training <- omic_methylation_training[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing methylation traning dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(omic_methylation_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_methylation_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_methylation_training <- rownames(omic_methylation_training)[idx]

nrow(omic_methylation_training)
ncol(omic_methylation_training)

## matching participants with data for training_clinical_numeric and omic_methylation_training

common_methylation_ids_training <- rownames(training_clinical_numeric)
common_methylation_ids_training <- intersect(common_methylation_ids_training, colnames(omic_methylation_training))

## construct a dataset with data on methylation (in transpose) and clinical outcome

dataset_methylation_training <- t(omic_methylation_training[predictors_omic_methylation_training,common_methylation_ids_training])

class(dataset_methylation_training)

# investigating the presence of missing data
d <- data.frame(dataset_methylation_training)
miss_var_summary(d)

# no missingness is observed in dataset_methylation_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_methylation_training <- training_clinical_numeric[common_methylation_ids_training,outcome.var]

outcome_methylation_training <- as.factor(outcome_methylation_training)
outcome_methylation_training <- ifelse(outcome_methylation_training==0,"No","Yes")


## standardize features (variance and scaling) 
dataset_methylation_training_var <- apply(dataset_methylation_training,2,var)
dataset_methylation_training <- dataset_methylation_training[,dataset_methylation_training_var > 2e-16]
dataset_methylation_training <- scale(dataset_methylation_training)

class(dataset_methylation_training)

## specify controls and train the models
model_list_methylation <- train.model(dataset_methylation_training, outcome_methylation_training)

########################################
# train_model_list
resamp <- resamples(model_list_methylation)

summary(resamp)

# plot the models
lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_2_6/fig_roc_plot_methylation.jpeg")

dotplot(resamp)

# rf model
dataset_methylation_training <- cbind(dataset_methylation_training, outcome_methylation_training)


train.predicts_methylation <- predict(
  model_list_methylation[["rf"]], 
  newdata = dataset_methylation_training)


cm <- confusionMatrix(train.predicts_methylation, as.factor(outcome_methylation_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.methylation <- data.frame(Model = "rf", 
                                     Accuracy = accuracy, 
                                     CI_95_percent = paste(ci_lower, "-" , ci_upper)
)
