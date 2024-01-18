#---
#title: "Model 2 testing - omic mutation"
#author: 'StduentNo : 2332635'
#---


src2 <- paste(wd, "/scripts/preprocessing/testing.R", sep = "")
source(src2)

dim(omic_mutation_testing)

## remove features with no variance or all missing values

feature.var <- apply(omic_mutation_testing,1,var,na.rm=T)
omic_mutation_testing <- omic_mutation_testing[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mutation testing dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(omic_mutation_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mutation_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mutation_testing <- rownames(omic_mutation_testing)[idx]

nrow(omic_mutation_testing)
ncol(omic_mutation_testing)

## matching participants with data for testing_clinical_numeric and omic_mutation_testing

common_mutation_ids_testing <- rownames(testing_clinical_numeric)
common_mutation_ids_testing <- intersect(common_mutation_ids_testing, colnames(omic_mutation_testing))

## construct a dataset with data on mutation (in transpose) and clinical outcome

dataset_mutation_testing <- t(omic_mutation_testing[predictors_omic_mutation_training,common_mutation_ids_testing])

class(dataset_mutation_testing)

# investigating the presence of missing data
d <- data.frame(dataset_mutation_testing)
miss_var_summary(d)

# no missingness is observed in dataset_mutation_clinical as well as in the output, pfi
# as no missingness, imputation is excluded


## outcome variable
outcome_mutation_testing <- testing_clinical_numeric[common_mutation_ids_testing,outcome.var]

outcome_mutation_testing <- as.factor(outcome_mutation_testing)
outcome_mutation_testing <- ifelse(outcome_mutation_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_mutation_testing_var <- apply(dataset_mutation_testing,2,var)
dataset_mutation_testing <- dataset_mutation_testing[,dataset_mutation_testing_var > 2e-16]
dataset_mutation_testing <- scale(dataset_mutation_testing)


#prediction
# rf model
prediction <- predict(
  model_list[["rf"]],
  newdata = dataset_mutation_testing)

confusionMatrix(prediction, as.factor(outcome_mutation_testing))
