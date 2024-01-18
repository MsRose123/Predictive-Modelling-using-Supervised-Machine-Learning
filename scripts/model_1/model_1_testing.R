#---
#title: "Model 1 - testing dataset - clinical features"
#author: 'StduentNo : 2332635'
#---

src <- paste(wd, "/scripts/preprocessingFunctions.R", sep = "")
source(src)

# explore testing_clinical dataset

#summary
head(testing_clinical)
summary(testing_clinical)

# distinct values
distinct_values(testing_clinical) 

# remove pfi.time and store in a new variable
testing_clinical_numeric <- select(testing_clinical, -pfi.time)

# converting each non-numeric variable to numeric 

# histology
# respecifiying NA as null values as those are stored as character in this feature
testing_clinical_numeric$histology <- ifelse(testing_clinical_numeric$histology == "NA", NA, testing_clinical_numeric$histology)
testing_clinical_numeric$histology <- as.numeric(as.factor(testing_clinical_numeric$histology))

# her2.status
testing_clinical_numeric$her2.status <- as.numeric(as.factor(testing_clinical_numeric$her2.status))

# ethnicity
testing_clinical_numeric$ethnicity <- ifelse(testing_clinical_numeric$ethnicity == "NA", NA, testing_clinical_numeric$ethnicity)
testing_clinical_numeric$ethnicity <- as.numeric(as.factor(testing_clinical_numeric$ethnicity))

# race
testing_clinical_numeric$race <- ifelse(testing_clinical_numeric$race == "NA", NA, testing_clinical_numeric$race)
testing_clinical_numeric$race <- as.numeric(as.factor(testing_clinical_numeric$race))

# stage
testing_clinical_numeric$stage <- ifelse(testing_clinical_numeric$stage == "NA", NA, testing_clinical_numeric$stage)
testing_clinical_numeric$stage <- as.numeric(as.factor(testing_clinical_numeric$stage))

# tnm.m.category
testing_clinical_numeric$tnm.m.category <- as.numeric(as.factor(testing_clinical_numeric$tnm.m.category))

# tnm.n.category
testing_clinical_numeric$tnm.n.category <- as.numeric(as.factor(testing_clinical_numeric$tnm.n.category))

# tnm.t.category
testing_clinical_numeric$tnm.t.category <- as.numeric(as.factor(testing_clinical_numeric$tnm.t.category))

# estrogen.receptor.status
testing_clinical_numeric$estrogen.receptor.status <- as.numeric(as.factor(testing_clinical_numeric$estrogen.receptor.status))

# progesterone.receptor.status
testing_clinical_numeric$progesterone.receptor.status <- as.numeric(as.factor(testing_clinical_numeric$progesterone.receptor.status))

# distinct values
distinct_values(testing_clinical_numeric)

# investigating the presence of missing data
#missing_data(training_clinical_numeric)
miss_var_summary(testing_clinical_numeric)

## outcome : progression-free interval
outcome_clinical_testing <- as.numeric(testing_clinical_numeric[,outcome.var])

# Remove the target variable before imputing ( pfi does not have missingness)
testing_clinical_numeric_1 <- select(testing_clinical_numeric, -pfi)

## explore the presence of outliers in the dataset using IQR
outliers <- presence.of.outliers(testing_clinical_numeric_1)

## impute missing values with the median value for the feature
testing_clinical_numeric_imputed <- impute.data(testing_clinical_numeric_1)

# investigating the presence of missing data after imputing
miss_var_summary((testing_clinical_numeric_imputed))

## standardize features
class(testing_clinical_numeric_imputed)

clinical_imp <- preProcess(testing_clinical_numeric_imputed, method=c('center', 'scale'))
testing_clinical_numeric_imputed <- predict(clinical_imp, testing_clinical_numeric_imputed)

# updated dataset with features selected by boruta feature selection method
testing_clinical_dataset <- testing_clinical_numeric_imputed[,clinical.vars]

nrow(testing_clinical_dataset)
ncol(testing_clinical_dataset)

# combine the data and convert to matrix
testing_clinical_dataset <- do.call(cbind, testing_clinical_dataset)

# run PCA using prcomp to remove the collinarity between the data 

pcfit <- prcomp(testing_clinical_dataset)

testing_clinical_dataset <- pcfit$x

## processing outcome variable
testing_clinical_outcome <- as.factor(outcome_clinical_testing)
testing_clinical_outcome <- ifelse(testing_clinical_outcome==0,"No","Yes")

# correcting the levels
testing_clinical_dataset_new <- as.data.frame(testing_clinical_dataset)                              # Duplicate test data set
training_clinical_dataset_new <- as.data.frame(training_clinical_dataset)

testing_clinical_dataset_new$PC1[which(!(testing_clinical_dataset_new$PC1 %in% unique(training_clinical_dataset_new$PC1)))] <- NA  # Replace new levels by NA

# convert to matrix
testing_clinical_dataset_new <- do.call(cbind, testing_clinical_dataset_new)

class(model_list)
class(training_clinical_dataset_new)
#model_list_new <- subset(model_list, select = -rf)

#model_list_new <- within(model_list, rm(rf))

## predict the model outcome

testing_clinical_dataset_new <- cbind(testing_clinical_dataset_new, testing_clinical_outcome)

model.predict.clinical <- predict(
  model_list[["rf"]], 
  newdata = testing_clinical_dataset_new)


# run confusion matrix for the model with the highest accuracy here
confusionMatrix(model.predict.clinical, testing_clinical_outcome) 



############## add graphs to represent data at each point
