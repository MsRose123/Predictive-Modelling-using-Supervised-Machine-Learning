#---
#title: "Preprocessing Functions"
#author: 'StduentNo : 2332635'
#---


#############################

## function : my.read.table
## for loading tab-delimited spreadsheets

library(data.table)
my.read.table <- function(filename, ...) {
  cat("reading", basename(filename), "... ")
  ## read in tab-delimited spreadsheet
  x <- fread(
    filename,
    header=T,
    stringsAsFactors=F,
    sep="\t",
    ...)
  ## remove any duplicate rows (identified by the first column)
  x <- x[match(unique(x[[1]]), x[[1]]),]
  ## make the first column the rownames of the data frame
  x <- data.frame(x,row.names=1,stringsAsFactors=F)
  cat(nrow(x), "x", ncol(x), "\n")
  x
}

#############################

## function : read.data
# function to read data in training and testing paths

read.data <- function(path.dir){
  
    ## list all files in the training dataset
    filenames <- list.files(path.dir, full.names=T)
    
    ## for this example we'll omit CNV
    filenames <- filenames[!grepl("cnv", filenames)]
    
    # files that does not include example, annotation, targets.
    dat.filenames <- filenames[-grep("(example|annotation|targets)", filenames)]
    
    ## load the data files into a list
    data0 <- lapply(dat.filenames, my.read.table)
    
    ## name the items of the list by the filename
    names(data0) <- sub(".txt", "", basename(dat.filenames))
    data0
}

#############################

## function : clinical.remove
# function to remove clinical from each dataset

#clinical.remove <- function(data0){    
    ## remove the clinical data from the list,
 #   clinical <- data0$clinical
  #  data0$clinical <- NULL
   # clinical

#}

#############################

## outcome if progression-free interval
outcome.var <- "pfi"

#############################

## function to obtain distinct value for each column in a data frame

distinct_values <- function(dataset){
  lapply(dataset, unique)
}

#############################

# investigating the presence of missing data
missing_data <- function(dataset) {
  sapply(dataset, function(x) sum(is.na(x)))
}

#############################

# explore the presence of outliers in the dataset using IQR

presence.of.outliers <- function(dataset){
  
  outliers_list <- list()
  
  for (var in names(dataset)) {
    IQR <- IQR(dataset[[var]], na.rm = TRUE)
    
    lower_bound <- quantile(dataset[[var]], 0.25, na.rm = TRUE) - 1.5 * IQR
    upper_bound <- quantile(dataset[[var]], 0.75, na.rm = TRUE) + 1.5 * IQR
    
    outliers <- dataset[[var]] < lower_bound | dataset[[var]] > upper_bound
    
    outliers_list[[var]] <- dataset[outliers, ]
    
  }
  outliers_list
}

#############################

## impute missing values with the median value for the feature
impute.data <- function(dataset){
    idx <- which(is.na(dataset),arr.ind=T)
    median.values <- apply(dataset,2,median,na.rm=T)
    dataset[idx] <- median.values[idx[,2]]
    dataset
}

#############################

## standardize features

standardise.features <- function(dataset){
  feature.var <- apply(dataset,2,var)
  dataset <- dataset[,feature.var > 2e-16]
  dataset <- scale(dataset)
  dataset
}

#############################

