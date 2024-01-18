#---
#title: "Model 8 training - omic gene and protein combined with clinical features"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

# protein and gene combined model
dataset_model_8_training <- combined_df_omics_training

## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna and protein omics combined with clinical features", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(dataset_model_8_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(dataset_model_8_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_model_8 <- rownames(dataset_model_8_training)[idx]

nrow(dataset_model_8_training)
ncol(dataset_model_8_training)

## matching participants with data for training_clinical_numeric and dataset_model_8_training

common_model_8_ids_training <- rownames(training_clinical_numeric)
common_model_8_ids_training <- intersect(common_model_8_ids_training, colnames(dataset_model_8_training))

## construct a dataset with data on dataset_model_8_training (in transpose)

dataset_combined_model_8_training <- t(dataset_model_8_training[predictors_model_8,common_model_8_ids_training])

class(dataset_combined_model_8_training)
nrow(dataset_combined_model_8_training)
ncol(dataset_combined_model_8_training)

# convert matrix to df
dataset_combined_model_8_training <- data.frame(dataset_combined_model_8_training)
dataset_combined_model_8_training

# adding columns from training_clinical_numeric as a df at the end of dataset_combined_model_8
dataset_combined_model_8_training$training_clinical_numeric <- training_clinical_numeric[common_model_8_ids_training,setdiff(colnames(training_clinical_numeric),outcome.var)]


## merge data for each data type into a single matrix
dataset_combined_model_8_training <- do.call(cbind, dataset_combined_model_8_training)

# investigating the presence of missing data
d <- data.frame(dataset_combined_model_8_training)
miss_var_summary(d)

## impute missing values with the median value for the feature
idx <- which(is.na(dataset_combined_model_8_training),arr.ind=T)
median.values <- apply(dataset_combined_model_8_training,2,median,na.rm=T)
dataset_combined_model_8_training[idx] <- median.values[idx[,2]]

miss_var_summary(dataset_combined_model_8_training)

## outcome variable
outcome_model_8_training <- training_clinical_numeric[common_model_8_ids_training,outcome.var]

outcome_model_8_training <- as.factor(outcome_model_8_training)

outcome_model_8_training <- ifelse(outcome_model_8_training==0,"No","Yes")

## standardize features (variance and scaling) 
dataset_combined_model_8_training_var <- apply(dataset_combined_model_8_training,2,var)
dataset_combined_model_8_training <- dataset_combined_model_8_training[,dataset_combined_model_8_training_var > 2e-16]
dataset_combined_model_8_training <- scale(dataset_combined_model_8_training)

## specify controls and train the models
model_list_combined_clinical_omics <- train.model(dataset_combined_model_8_training, outcome_model_8_training)

########################################
# train_model_list
resamp <- resamples(model_list_combined_clinical_omics)

summary(resamp)

# plot the models

lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_8/fig_roc_plot_combined_clinical_omics.jpeg")

dotplot(resamp)

# rf model
dataset_combined_model_8_training <- cbind(dataset_combined_model_8_training, outcome_model_8_training)


train.predicts_combined_clinical_omics <- predict(
  model_list_combined_clinical_omics[["rf"]], 
  newdata = dataset_combined_model_8_training)


cm <- confusionMatrix(train.predicts_combined_clinical_omics, as.factor(outcome_model_8_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.model_8 <- data.frame(Model = "rf", 
                                     Accuracy = accuracy, 
                                     CI_95_percent = paste(ci_lower, "-" , ci_upper)
)
