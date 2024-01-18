#---
#title: "Model 2 testing - omic protein"
#author: 'StduentNo : 2332635'
#---


src2 <- paste(wd, "/scripts/preprocessing/testing.R", sep = "")
source(src2)

dim(omic_protein_testing)

## remove features with no variance or all missing values

feature.var <- apply(omic_protein_testing,1,var,na.rm=T)
omic_protein_testing <- omic_protein_testing[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing protein testing dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(omic_protein_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_protein_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_protein_testing <- rownames(omic_protein_testing)[idx]

nrow(omic_protein_testing)
ncol(omic_protein_testing)

## matching participants with data for testing_clinical_numeric and omic_protein_testing

common_protein_ids_testing <- rownames(testing_clinical_numeric)
common_protein_ids_testing <- intersect(common_protein_ids_testing, colnames(omic_protein_testing))

## construct a dataset with data on protein (in transpose) and clinical outcome

dataset_protein_testing <- t(omic_protein_testing[predictors_omic_protein_training,common_protein_ids_testing])

class(dataset_protein_testing)

# investigating the presence of missing data
d <- data.frame(dataset_protein_testing)
miss_var_summary(d)

# no missingness is observed in dataset_protein_clinical as well as in the output, pfi
# as no missingness, imputation is excluded


## outcome variable
outcome_protein_testing <- testing_clinical_numeric[common_protein_ids_testing,outcome.var]

outcome_protein_testing <- as.factor(outcome_protein_testing)
outcome_protein_testing <- ifelse(outcome_protein_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_protein_testing_var <- apply(dataset_protein_testing,2,var)
dataset_protein_testing <- dataset_protein_testing[,dataset_protein_testing_var > 2e-16]
dataset_protein_testing <- scale(dataset_protein_testing)


#prediction
# rf model
prediction <- predict(
  model_list[["rf"]],
  newdata = dataset_protein_testing)

confusionMatrix(prediction, as.factor(outcome_protein_testing))
