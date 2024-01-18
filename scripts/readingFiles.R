#---
#title: "Reading files and loading data"
#author: 'StduentNo : 2332635'
#---

src <- paste(wd, "/scripts/preprocessing/preprocessingFunctions.R", sep = "")
source(src)

############# Training data

## list all files in the training dataset
filenames <- list.files(training.dir, full.names=T)

## for this example we'll omit CNV
filenames <- filenames[!grepl("cnv", filenames)]

dat.filenames <- filenames[-grep("(example|annotation|targets)", filenames)]

## load the data files into a list
training_data0 <- lapply(dat.filenames, my.read.table)

## name the items of the list by the filename
names(training_data0) <- sub(".txt", "", basename(dat.filenames))


## storing clinical data after removing from the list

training_clinical <- training_data0$clinical
training_data0$clinical <- NULL

# saving protein annotation files for future use
file_annotation_protein <- filenames[grepl("protein-annotation", filenames)]

############ Testing data

## list all files in the testing dataset
filenames <- list.files(testing.dir, full.names=T)

## omitting cnv
filenames <- filenames[!grepl("cnv", filenames)]

## load the data files into a list
testing_data0 <- lapply(filenames, my.read.table)

## name the items of the list by the filename
names(testing_data0) <- sub(".txt", "", basename(filenames))

## storing clinical data after removing from the list

testing_clinical <- testing_data0$clinical
testing_data0$clinical <- NULL

