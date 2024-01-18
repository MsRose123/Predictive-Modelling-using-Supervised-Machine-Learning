#---
#title: "Model 2 training - omic mirna"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

dim(omic_mirna_training)

## remove features with no variance or all missing values

feature.var <- apply(omic_mirna_training,1,var,na.rm=T)
omic_mirna_training <- omic_mirna_training[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mirna traning dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(omic_mirna_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mirna_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mirna_training <- rownames(omic_mirna_training)[idx]

nrow(omic_mirna_training)
ncol(omic_mirna_training)

## matching participants with data for training_clinical_numeric and omic_mirna_training

common_mirna_ids_training <- rownames(training_clinical_numeric)
common_mirna_ids_training <- intersect(common_mirna_ids_training, colnames(omic_mirna_training))

## construct a dataset with data on mirna (in transpose) and clinical outcome

dataset_mirna_training <- t(omic_mirna_training[predictors_omic_mirna_training,common_mirna_ids_training])

class(dataset_mirna_training)

# investigating the presence of missing data
d <- data.frame(dataset_mirna_training)
miss_var_summary(d)

# no missingness is observed in dataset_mirna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_mirna_training <- training_clinical_numeric[common_mirna_ids_training,outcome.var]

outcome_mirna_training <- as.factor(outcome_mirna_training)
outcome_mirna_training <- ifelse(outcome_mirna_training==0,"No","Yes")


## standardize features (variance and scaling) 
dataset_mirna_training_var <- apply(dataset_mirna_training,2,var)
dataset_mirna_training <- dataset_mirna_training[,dataset_mirna_training_var > 2e-16]
dataset_mirna_training <- scale(dataset_mirna_training)

class(dataset_mirna_training)

## specify controls and train the models
model_list_mirna <- train.model(dataset_mirna_training, outcome_mirna_training)

# train_model_list
resamp <- resamples(model_list_mirna)

summary(resamp)

# plot and save the models

lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_2_6/fig_roc_plot_mirna.jpeg")


dotplot(resamp)

# rf model
dataset_mirna_training <- cbind(dataset_mirna_training, outcome_mirna_training)

prediction <- predict(
  model_list_mirna[["rf"]],
  newdata = dataset_mirna_training)


cm <- confusionMatrix(prediction, as.factor(outcome_mirna_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.mirna <- data.frame(Model = "rf", 
                                      Accuracy = accuracy, 
                                      CI_95_percent = paste(ci_lower, "-" , ci_upper)
)


