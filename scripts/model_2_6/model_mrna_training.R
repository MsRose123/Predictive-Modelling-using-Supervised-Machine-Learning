#---
#title: "Model 2 training - omic mrna"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

dim(omic_mrna_training)

## remove features with no variance or all missing values

feature.var <- apply(omic_mrna_training,1,var,na.rm=T)
omic_mrna_training <- omic_mrna_training[which(feature.var > 2e-16),]


## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna traning dataset", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(omic_mrna_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(omic_mrna_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omic_mrna_training <- rownames(omic_mrna_training)[idx]

nrow(omic_mrna_training)
ncol(omic_mrna_training)

## matching participants with data for training_clinical_numeric and omic_mrna_training

common_mrna_ids_training <- rownames(training_clinical_numeric)
common_mrna_ids_training <- intersect(common_mrna_ids_training, colnames(omic_mrna_training))

## construct a dataset with data on mrna (in transpose) and clinical outcome

dataset_mrna_training <- t(omic_mrna_training[predictors_omic_mrna_training,common_mrna_ids_training])

class(dataset_mrna_training)

# investigating the presence of missing data
d <- data.frame(dataset_mrna_training)
miss_var_summary(d)

# no missingness is observed in dataset_mrna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_mrna_training <- training_clinical_numeric[common_mrna_ids_training,outcome.var]

outcome_mrna_training <- as.factor(outcome_mrna_training)
outcome_mrna_training <- ifelse(outcome_mrna_training==0,"No","Yes")


## standardize features (variance and scaling) 
dataset_mrna_training_var <- apply(dataset_mrna_training,2,var)
dataset_mrna_training <- dataset_mrna_training[,dataset_mrna_training_var > 2e-16]
dataset_mrna_training <- scale(dataset_mrna_training)

class(dataset_mrna_training)

## specify controls and train the models
model_list_mrna <- train.model(dataset_mrna_training, outcome_mrna_training)

########################################
# train_model_list
resamp <- resamples(model_list_mrna)

summary(resamp)

# plot the models

lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_2_6/fig_roc_plot_mrna.jpeg")

dotplot(resamp)

# rf model
dataset_mrna_training <- cbind(dataset_mrna_training, outcome_mrna_training)


train.predicts_mrna <- predict(
  model_list_mrna[["rf"]], 
  newdata = dataset_mrna_training)


cm <- confusionMatrix(train.predicts_mrna, as.factor(outcome_mrna_training))


# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.mrna <- data.frame(Model = "rf", 
                                   Accuracy = accuracy, 
                                   CI_95_percent = paste(ci_lower, "-" , ci_upper)
)

