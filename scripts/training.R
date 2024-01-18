#---
#title: "Traning script"
#author: 'StduentNo : 2332635'
#---

getwd()


# number of folds, k for cross-validation
#numFolds <- createFolds(outcome, k = 10)
#numFolds

# Verify that the folds maintain the proportion of yes/no results
#sapply(numFolds, function(i) {
#  table(outcome[i])/length(i)
#})


# training controls - common parameters for training diferent models 
myControl <- trainControl(
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = FALSE,
  savePredictions = TRUE,
  method="cv", 
  number=10
)

train.model <- function(train.dataset, train.outcome){
  
  # adding outcome back to the dataset and train the models
  train.dataset <- cbind(train.dataset, train.outcome)
  
  # elastic net model
  print("training elastic net model... ")
  glm_model <- train(train.outcome ~ .,
                     train.dataset,
                     method = "glmnet",
                     family = "binomial",
                     metric = "ROC",
                     tuneGrid = expand.grid(
                       alpha = 0:1,
                       lambda = 0:10/10),
                     trControl = myControl
  )
  # plot and save
  
  # random forest
  print("training random forest model... ")
  rf_model <- train(train.outcome ~ .,
                    train.dataset,
                    metric = "ROC",
                    method = "ranger",
                    tuneGrid = expand.grid(
                      mtry = c(2, 5, 10, 19),
                      splitrule = c("gini", "extratrees"),
                      min.node.size = 1),
                    trControl = myControl)
  # plot and save
  
  # knn model
  print("training knn model... ")
  knn_model <- train(train.outcome ~ .,
                     train.dataset,
                     metric = "ROC",
                     method = "knn",
                     tuneLength = 20,
                     trControl = myControl)
  # plot and save
  
  # svm
  print("training svm model... ")
  svm_model <- train(train.outcome ~ .,
                     train.dataset,
                     metric = "ROC",
                     method = "svmRadial",
                     tuneLength = 10,
                     trControl = myControl,
                     maxit = 1000
                     )
  # plot and save
  
  # naive bayes
  print("training naive bayes model... ")
  nb_model <- train(train.outcome ~ .,
                    train.dataset,
                    metric = "ROC",
                    method = "naive_bayes",
                    trControl = myControl)
  # plot and save
  
  # comparing model fitting
  print("comparing model fitting... ")
  model_list <- list(glmmet = glm_model,
                     rf = rf_model,
                     knn = knn_model,
                     svm = svm_model,
                     nb = nb_model)
  
  model_list
}
