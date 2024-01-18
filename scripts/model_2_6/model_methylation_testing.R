#---
#title: "Model 2 testing - omic methylation"
#author: 'StduentNo : 2332635'
#---


src2 <- paste(wd, "/scripts/preprocessing/testing.R", sep = "")
source(src2)

dim(omic_methylation_testing)

## remove features with no variance or all missing values

feature.var <- apply(omic_methylation_testing,1,var,na.rm=T)
omic_methylation_testing <- omic_methylation_testing[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing methylation testing dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(omic_methylation_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_methylation_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_methylation_testing <- rownames(omic_methylation_testing)[idx]

nrow(omic_methylation_testing)
ncol(omic_methylation_testing)

## matching participants with data for testing_clinical_numeric and omic_methylation_testing

common_methylation_ids_testing <- rownames(testing_clinical_numeric)
common_methylation_ids_testing <- intersect(common_methylation_ids_testing, colnames(omic_methylation_testing))

## construct a dataset with data on methylation (in transpose) and clinical outcome

dataset_methylation_testing <- t(omic_methylation_testing[predictors_omic_methylation_training,common_methylation_ids_testing])

class(dataset_methylation_testing)

# investigating the presence of missing data
d <- data.frame(dataset_methylation_testing)
miss_var_summary(d)

# no missingness is observed in dataset_methylation_clinical as well as in the output, pfi
# as no missingness, imputation is excluded


## outcome variable
outcome_methylation_testing <- testing_clinical_numeric[common_methylation_ids_testing,outcome.var]

outcome_methylation_testing <- as.factor(outcome_methylation_testing)
outcome_methylation_testing <- ifelse(outcome_methylation_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_methylation_testing_var <- apply(dataset_methylation_testing,2,var)
dataset_methylation_testing <- dataset_methylation_testing[,dataset_methylation_testing_var > 2e-16]
dataset_methylation_testing <- scale(dataset_methylation_testing)


#prediction
# rf model
# train_model_list
resamp <- resamples(model_list_methylation)

summary(resamp)

# plot the models
lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_2_6/fig_roc_plot_methylation.jpeg")

dotplot(resamp)

# rf model

train.predicts_methylation <- predict(
  model_list_methylation[["rf"]], 
  newdata = dataset_methylation_testing)

confusionMatrix(train.predicts_methylation, as.factor(outcome_methylation_testing))

