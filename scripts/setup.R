#---
#title: "setup"
#author: 'StduentNo : 2332635'
#---

# clear environment
rm(list=ls())

# system date
Sys.Date()

# set working directory
wd <- getwd()
setwd(wd)


# Load necessary packages and libraries

install.packages("dplyr") #data manipulation
install.packages("naniar") # missing data
install.packages("magrittr")
install.packages("precrec")
install.packages("visdat")
install.packages("data.table")
install.packages("mlr3verse")
install.packages("FactoMineR")
install.packages("missMDA")
install.packages("ggrepel")
install.packages("limma")
install.packages("BiocManager")
install.packages("glmnet")
install.packages("naivebayes")
install.packages("caret")
install.packages("e1071")
install.packages("Boruta") #feature selection
install.packages("lattice")

library(dplyr)
library(naniar)
library(gtsummary)
library(magrittr) ## for '%>%' operator
library(precrec) ## for mlr3 plots
library(visdat)
library(data.table)
library(mlr3verse)
library(FactoMineR)
library(missMDA)
library(ggrepel)
library(limma)
library(glmnet)
library(naivebayes)
library(caret)
library(e1071)
library(Boruta)
library(lattice)


# load training and testing dataset
training.dir <- paste(wd, "/data/trainingDataset/rawData", sep = "")
testing.dir <- paste(wd, "/data/testingDataset/rawData", sep = "")


# model 1 - clinical features 
# model 2 - mirna features
# model 3 - mrna features
# model 4 - mutation
# model 5 - protein
# model 6 - methylation
# model 7 - features from two omics - mrna and protein
# model 8 - features from clinical and two omic (mrna and protein)

# record of computational environment and save it in a file
#sessionInfo()
