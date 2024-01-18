#---
#title: "Model 2 testing - omic mirna"
#author: 'StduentNo : 2332635'
#---


src2 <- paste(wd, "/scripts/preprocessing/testing.R", sep = "")
source(src2)

dim(omic_mirna_testing)

## remove features with no variance or all missing values

feature.var <- apply(omic_mirna_testing,1,var,na.rm=T)
omic_mirna_testing <- omic_mirna_testing[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mirna testing dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(omic_mirna_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mirna_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mirna_testing <- rownames(omic_mirna_testing)[idx]

nrow(omic_mirna_testing)
ncol(omic_mirna_testing)

## matching participants with data for testing_clinical_numeric and omic_mirna_testing

common_mirna_ids_testing <- rownames(testing_clinical_numeric)
common_mirna_ids_testing <- intersect(common_mirna_ids_testing, colnames(omic_mirna_testing))

## construct a dataset with data on mirna (in transpose) and clinical outcome

dataset_mirna_testing <- t(omic_mirna_testing[predictors_omic_mirna_training,common_mirna_ids_testing])

class(dataset_mirna_testing)

# investigating the presence of missing data
d <- data.frame(dataset_mirna_testing)
miss_var_summary(d)

# no missingness is observed in dataset_mirna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded


## outcome variable
outcome_mirna_testing <- testing_clinical_numeric[common_mirna_ids_testing,outcome.var]

outcome_mirna_testing <- as.factor(outcome_mirna_testing)
outcome_mirna_testing <- ifelse(outcome_mirna_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_mirna_testing_var <- apply(dataset_mirna_testing,2,var)
dataset_mirna_testing <- dataset_mirna_testing[,dataset_mirna_testing_var > 2e-16]
dataset_mirna_testing <- scale(dataset_mirna_testing)


#prediction
# rf model
prediction <- predict(
  model_list[["rf"]],
  newdata = dataset_mirna_testing)

confusionMatrix(prediction, as.factor(outcome_mirna_testing))
