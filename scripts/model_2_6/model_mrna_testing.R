#---
#title: "Model 2 testing - omic mrna"
#author: 'StduentNo : 2332635'
#---


src2 <- paste(wd, "/scripts/preprocessing/testing.R", sep = "")
source(src2)

dim(omic_mrna_testing)

## remove features with no variance or all missing values

feature.var <- apply(omic_mrna_testing,1,var,na.rm=T)
omic_mrna_testing <- omic_mrna_testing[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna testing dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(omic_mrna_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mrna_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mrna_testing <- rownames(omic_mrna_testing)[idx]

nrow(omic_mrna_testing)
ncol(omic_mrna_testing)

## matching participants with data for testing_clinical_numeric and omic_mrna_testing

common_mrna_ids_testing <- rownames(testing_clinical_numeric)
common_mrna_ids_testing <- intersect(common_mrna_ids_testing, colnames(omic_mrna_testing))

## construct a dataset with data on mrna (in transpose) and clinical outcome

dataset_mrna_testing <- t(omic_mrna_testing[predictors_omic_mrna_training,common_mrna_ids_testing])

class(dataset_mrna_testing)

# investigating the presence of missing data
d <- data.frame(dataset_mrna_testing)
miss_var_summary(d)

# no missingness is observed in dataset_mrna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded


## outcome variable
outcome_mrna_testing <- testing_clinical_numeric[common_mrna_ids_testing,outcome.var]

outcome_mrna_testing <- as.factor(outcome_mrna_testing)
outcome_mrna_testing <- ifelse(outcome_mrna_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_mrna_testing_var <- apply(dataset_mrna_testing,2,var)
dataset_mrna_testing <- dataset_mrna_testing[,dataset_mrna_testing_var > 2e-16]
dataset_mrna_testing <- scale(dataset_mrna_testing)


#prediction
# rf model
prediction <- predict(
  model_list[["rf"]],
  newdata = dataset_mrna_testing)

confusionMatrix(prediction, as.factor(outcome_mrna_testing))
