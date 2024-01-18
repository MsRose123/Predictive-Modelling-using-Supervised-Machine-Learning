#---
#title: "Model 2 training - omic protein"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

dim(omic_protein_training)

## remove features with no variance or all missing values

feature.var <- apply(omic_protein_training,1,var,na.rm=T)
omic_protein_training <- omic_protein_training[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing protein traning dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(omic_protein_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_protein_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_protein_training <- rownames(omic_protein_training)[idx]

nrow(omic_protein_training)
ncol(omic_protein_training)

## matching participants with data for training_clinical_numeric and omic_protein_training

common_protein_ids_training <- rownames(training_clinical_numeric)
common_protein_ids_training <- intersect(common_protein_ids_training, colnames(omic_protein_training))

## construct a dataset with data on protein (in transpose) and clinical outcome

dataset_protein_training <- t(omic_protein_training[predictors_omic_protein_training,common_protein_ids_training])

class(dataset_protein_training)

# investigating the presence of missing data
d <- data.frame(dataset_protein_training)
miss_var_summary(d)

# no missingness is observed in dataset_protein_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_protein_training <- training_clinical_numeric[common_protein_ids_training,outcome.var]

outcome_protein_training <- as.factor(outcome_protein_training)
outcome_protein_training <- ifelse(outcome_protein_training==0,"No","Yes")


## standardize features (variance and scaling) 
dataset_protein_training_var <- apply(dataset_protein_training,2,var)
dataset_protein_training <- dataset_protein_training[,dataset_protein_training_var > 2e-16]
dataset_protein_training <- scale(dataset_protein_training)

class(dataset_protein_training)

## specify controls and train the models
model_list_protein <- train.model(dataset_protein_training, outcome_protein_training)

########################################
# train_model_list
resamp <- resamples(model_list_protein)

summary(resamp)

# plot the models

jpeg(file = "results/model_2_6/fig_roc_plot_protein.jpeg")
lattice::bwplot(resamp, metric = "ROC")

dotplot(resamp)

# rf model
dataset_protein_training <- cbind(dataset_protein_training, outcome_protein_training)


train.predicts_protein <- predict(
  model_list_protein[["rf"]], 
  newdata = dataset_protein_training)


cm <- confusionMatrix(train.predicts_protein, as.factor(outcome_protein_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.protein <- data.frame(Model = "rf", 
                                      Accuracy = accuracy, 
                                      CI_95_percent = paste(ci_lower, "-" , ci_upper)
)

