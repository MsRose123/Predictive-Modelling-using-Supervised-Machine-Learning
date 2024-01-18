#---
#title: "Model 8 testing - omic gene and protein combined with clinical features"
#author: 'StduentNo : 2332635'
#---

# protein and gene combined model
dataset_model_8_testing <- combined_df_omics_testing

## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna and protein omics combined with clinical features", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(dataset_model_8_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(dataset_model_8_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_model_8 <- rownames(dataset_model_8_testing)[idx]

nrow(dataset_model_8_testing)
ncol(dataset_model_8_testing)

## matching participants with data for testing_clinical_numeric and dataset_model_8_testing

common_model_8_ids_testing <- rownames(testing_clinical_numeric)
common_model_8_ids_testing <- intersect(common_model_8_ids_testing, colnames(dataset_model_8_testing))

## construct a dataset with data on dataset_model_8_testing (in transpose)

dataset_combined_model_8_testing <- t(dataset_model_8_testing[predictors_model_8,common_model_8_ids_testing])

class(dataset_combined_model_8_testing)
nrow(dataset_combined_model_8_testing)
ncol(dataset_combined_model_8_testing)

# convert matrix to df
dataset_combined_model_8_testing <- data.frame(dataset_combined_model_8_testing)
dataset_combined_model_8_testing

# adding columns from testing_clinical_numeric as a df at the end of dataset_combined_model_8
dataset_combined_model_8_testing$testing_clinical_numeric <- testing_clinical_numeric[common_model_8_ids_testing,setdiff(colnames(testing_clinical_numeric),outcome.var)]


## merge data for each data type into a single matrix
dataset_combined_model_8_testing <- do.call(cbind, dataset_combined_model_8_testing)

# investigating the presence of missing data
d <- data.frame(dataset_combined_model_8_testing)
miss_var_summary(d)

## impute missing values with the median value for the feature
idx <- which(is.na(dataset_combined_model_8_testing),arr.ind=T)
median.values <- apply(dataset_combined_model_8_testing,2,median,na.rm=T)
dataset_combined_model_8_testing[idx] <- median.values[idx[,2]]

miss_var_summary(dataset_combined_model_8_testing)

## outcome variable
outcome_model_8_testing <- testing_clinical_numeric[common_model_8_ids_testing,outcome.var]

outcome_model_8_testing <- as.factor(outcome_model_8_testing)

outcome_model_8_testing <- ifelse(outcome_model_8_testing==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_combined_model_8_testing_var <- apply(dataset_combined_model_8_testing,2,var)
dataset_combined_model_8_testing <- dataset_combined_model_8_testing[,dataset_combined_model_8_testing_var > 2e-16]
dataset_combined_model_8_testing <- scale(dataset_combined_model_8_testing)


# rf model
dataset_combined_model_8_testing <- cbind(dataset_combined_model_8_testing, outcome_model_8_testing)


train.predicts_combined_clinical_omics <- predict(
  model_list_combined_omics[["rf"]], 
  newdata = dataset_combined_model_8_testing)

confusionMatrix(train.predicts_combined_clinical_omics, as.factor(outcome_model_8_testing))


