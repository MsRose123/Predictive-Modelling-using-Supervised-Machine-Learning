#---
#title: "Model 2 training - omic mutation"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

dim(omic_mutation_training)

## remove features with no variance or all missing values

feature.var <- apply(omic_mutation_training,1,var,na.rm=T)
omic_mutation_training <- omic_mutation_training[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mutation traning dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(omic_mutation_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mutation_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mutation_training <- rownames(omic_mutation_training)[idx]

nrow(omic_mutation_training)
ncol(omic_mutation_training)

## matching participants with data for training_clinical_numeric and omic_mutation_training

common_mutation_ids_training <- rownames(training_clinical_numeric)
common_mutation_ids_training <- intersect(common_mutation_ids_training, colnames(omic_mutation_training))

## construct a dataset with data on mutation (in transpose) and clinical outcome

dataset_mutation_training <- t(omic_mutation_training[predictors_omic_mutation_training,common_mutation_ids_training])

class(dataset_mutation_training)

# investigating the presence of missing data
d <- data.frame(dataset_mutation_training)
miss_var_summary(d)

# no missingness is observed in dataset_mutation_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_mutation_training <- training_clinical_numeric[common_mutation_ids_training,outcome.var]

outcome_mutation_training <- as.factor(outcome_mutation_training)
outcome_mutation_training <- ifelse(outcome_mutation_training==0,"No","Yes")


## standardize features (variance and scaling) 
dataset_mutation_training_var <- apply(dataset_mutation_training,2,var)
dataset_mutation_training <- dataset_mutation_training[,dataset_mutation_training_var > 2e-16]
dataset_mutation_training <- scale(dataset_mutation_training)

class(dataset_mutation_training)

## specify controls and train the models
model_list_mutation <- train.model(dataset_mutation_training, outcome_mutation_training)

########################################
# train_model_list
resamp <- resamples(model_list_mutation)

summary(resamp)

# plot the models

jpeg(file = "results/model_2_6/fig_roc_plot_mutation.jpeg")
lattice::bwplot(resamp, metric = "ROC")

dotplot(resamp)

# rf model
dataset_mutation_training <- cbind(dataset_mutation_training, outcome_mutation_training)


train.predicts_mutation <- predict(
  model_list_mutation[["rf"]], 
  newdata = dataset_mutation_training)


cm <- confusionMatrix(train.predicts_mutation, as.factor(outcome_mutation_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.mutation <- data.frame(Model = "rf", 
                                  Accuracy = accuracy, 
                                  CI_95_percent = paste(ci_lower, "-" , ci_upper)
)
