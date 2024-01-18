#---
#title: "Model 1 - training dataset - clinical features"
#author: 'StduentNo : 2332635'
#---

src <- paste(wd, "/scripts/preprocessingFunctions.R", sep = "")
source(src)

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

# explore training_clinical dataset

#summary
head(training_clinical)
summary(training_clinical)

# distinct values
distinct_values(training_clinical) 

# remove pfi.time and store in a new variable
training_clinical_numeric <- select(training_clinical, -pfi.time)

# converting each non-numeric variable to numeric 

# histology
# respecifiying NA as null values as those are stored as character in this feature
training_clinical_numeric$histology <- ifelse(training_clinical_numeric$histology == "NA", NA, training_clinical_numeric$histology)
training_clinical_numeric$histology <- as.numeric(as.factor(training_clinical_numeric$histology))

# her2.status
training_clinical_numeric$her2.status <- as.numeric(as.factor(training_clinical_numeric$her2.status))

# ethnicity
training_clinical_numeric$ethnicity <- ifelse(training_clinical_numeric$ethnicity == "NA", NA, training_clinical_numeric$ethnicity)
training_clinical_numeric$ethnicity <- as.numeric(as.factor(training_clinical_numeric$ethnicity))

# race
training_clinical_numeric$race <- ifelse(training_clinical_numeric$race == "NA", NA, training_clinical_numeric$race)
training_clinical_numeric$race <- as.numeric(as.factor(training_clinical_numeric$race))

# stage
training_clinical_numeric$stage <- ifelse(training_clinical_numeric$stage == "NA", NA, training_clinical_numeric$stage)
training_clinical_numeric$stage <- as.numeric(as.factor(training_clinical_numeric$stage))

# tnm.m.category
training_clinical_numeric$tnm.m.category <- as.numeric(as.factor(training_clinical_numeric$tnm.m.category))

# tnm.n.category
training_clinical_numeric$tnm.n.category <- as.numeric(as.factor(training_clinical_numeric$tnm.n.category))

# tnm.t.category
training_clinical_numeric$tnm.t.category <- as.numeric(as.factor(training_clinical_numeric$tnm.t.category))

# estrogen.receptor.status
training_clinical_numeric$estrogen.receptor.status <- as.numeric(as.factor(training_clinical_numeric$estrogen.receptor.status))

# progesterone.receptor.status
training_clinical_numeric$progesterone.receptor.status <- as.numeric(as.factor(training_clinical_numeric$progesterone.receptor.status))

# distinct values
distinct_values(training_clinical_numeric)

# investigating the presence of missing data
#missing_data(training_clinical_numeric)
miss_var_summary(training_clinical_numeric)

## outcome : progression-free interval
outcome_clinical_training <- as.numeric(training_clinical_numeric[,outcome.var])

# number of 0's and 1's in the target variable
table(outcome)

# proportion of 0's and 1's in the target variable
table(outcome)/length(outcome)

# Remove the target variable before imputing ( pfi does not have missingness)
training_clinical_numeric_1 <- select(training_clinical_numeric, -pfi)

class(training_clinical_numeric_1)
## explore the presence of outliers in the dataset using IQR
outliers <- presence.of.outliers(training_clinical_numeric_1)

## impute missing values with the median value for the feature
training_clinical_numeric_imputed <- impute.data(training_clinical_numeric)

# investigating the presence of missing data after imputing
miss_var_summary((training_clinical_numeric_imputed))


## standardize features
#training_clinical_numeric_imputed <- standardise.features(training_clinical_numeric_imputed)

class(training_clinical_numeric_imputed)

clinical_imp <- preProcess(training_clinical_numeric_imputed, method=c('center', 'scale'))
training_clinical_numeric_imputed <- predict(clinical_imp, training_clinical_numeric_imputed)

############### feature selection

# Ensure the 'pfi' variable is a factor
training_clinical_numeric_imputed$pfi <- as.factor(outcome_clinical_training)

# Apply the Boruta feature selection method
boruta_output <- Boruta(pfi ~ ., data = training_clinical_numeric_imputed, doTrace = 0)
boruta_output

# Get the final decision about feature importance
clinical.vars <- getSelectedAttributes(boruta_output, withTentative = F)
clinical.vars

# Plot the boruta_output
plot(boruta_output, xlab = "", xaxt = "n")
############### 

# updated dataset with relevant features
training_clinical_dataset <- training_clinical_numeric_imputed[,clinical.vars]

class(training_clinical_dataset)

# convert to matrix
training_clinical_dataset <- do.call(cbind, training_clinical_dataset)

# run PCA using prcomp to remove the collinarity between the data 

pcfit <- prcomp(training_clinical_dataset)

training_clinical_dataset <- pcfit$x

## processing outcome variable
training_clinical_outcome <- as.factor(outcome_clinical_training)
training_clinical_outcome <- ifelse(training_clinical_outcome==0,"No","Yes")

## specify controls and train the models
model_list_clinical <- train.model(training_clinical_dataset, training_clinical_outcome)

class(training_clinical_dataset)
class(training_clinical_outcome)


# train_model_list
resamp <- resamples(model_list_clinical)

summary(resamp)

# plot the models
lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_1/fig_roc_plot_clinical.jpeg")


dotplot(resamp)

# rf model
training_clinical_dataset <- cbind(training_clinical_dataset, training_clinical_outcome)


train.predicts_clinical <- predict(
  model_list_clinical[["rf"]], 
  newdata = training_clinical_dataset)

cm <- confusionMatrix(train.predicts_clinical, as.factor(training_clinical_outcome))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.clinical <- data.frame(Model = "rf", 
                      Accuracy = accuracy, 
                      CI_95_percent = paste(ci_lower, "-" , ci_upper)
)

