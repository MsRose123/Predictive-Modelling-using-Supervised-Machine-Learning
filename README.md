# Predicting Breast Cancer Progression based on Multi-Omic and Clinical Data using Machine Learning Approach

The research aims to predict the progression of breast cancer using clinical information and multi-omic data. To compare different Machine Learning models built on combinations of different omics and clinical features (8 models built and trained on 5 ML approaches, KNN, SVM, Random Forest, naïve bayes, Logistic Regression, and Elastic Net algorithms). 

The Cancer Genome Atlas Program (TCGAP) database will be used. The research will employ median imputations to handle missing data and feature selection for each model will be done using the Boruta feature selection method. The models will be trained using the training dataset, and the predictive power of each model will be obtained by testing using the testing dataset. Each machine learning approach will be compared on the basis of prediction accuracy to find the best approach.

The application can be run using the shell script executable.sh
Run this on command line :

        $./executable.sh

Before running the pipeline add thhe testing and training datasets to the folder. The paths for this are as explained below.


Directories :
-------------

1. training.dir - pointing to the training dataset.
2. testing.dir - pointing to the testing dataset.

Data paths :
------------

1. training data - /data/trainingDataset/rawData
2. testing data - /data/testingDataset/rawData

Scripts :
---------
The scripts required for preprocessing, training and testing are stored in appropriate files and folders in the script folder

Results :
---------

The plots obtained and used in the report are stored in appropriate files and folders in the results folder

Directory structure :
---------------------

```
MLDM_assessment_2332635
├─ .DS_Store
├─ .RData
├─ .Rhistory
├─ .git
├─ .gitignore
├─ .snakemake
├─ MLDM_assessment_report_2332635.pdf
├─ README.md
├─ Snakefile
├─ data
│  ├─ .DS_Store
│  ├─ testingDataset
│  │  ├─ .DS_Store
│  │  └─ rawData
│  │     ├─ clinical.txt
│  │     ├─ cnv.txt
│  │     ├─ methylation.txt
│  │     ├─ mirna.txt
│  │     ├─ mrna.txt
│  │     ├─ mutation.txt
│  │     └─ protein.txt
│  └─ trainingDataset
│     ├─ .DS_Store
│     └─ rawData
│        ├─ .DS_Store
│        ├─ clinical.txt
│        ├─ cnv.txt
│        ├─ example.r
│        ├─ gene-annotation.txt
│        ├─ methylation-annotation.txt
│        ├─ methylation.txt
│        ├─ mirna-targets.txt
│        ├─ mirna.txt
│        ├─ mrna.txt
│        ├─ mutation.txt
│        ├─ protein-annotation.txt
│        └─ protein.txt
├─ environment.yml
├─ executable.sh
├─ results
│  ├─ .DS_Store
│  ├─ model_1
│  │  ├─ .DS_Store
│  │  └─ fig_roc_plot_clinical.jpeg
│  ├─ model_2_6
│  │  ├─ .DS_Store
│  │  ├─ fig_roc_plot_methylation.jpeg
│  │  ├─ fig_roc_plot_mirna.jpeg
│  │  ├─ fig_roc_plot_mrna.jpeg
│  │  ├─ fig_roc_plot_mutation.jpeg
│  │  └─ fig_roc_plot_protein.jpeg
│  ├─ model_7
│  │  └─ fig_roc_plot_combined_omics.jpeg
│  └─ model_8
│     ├─ .DS_Store
│     └─ fig_roc_plot_combined_clinical_omics.jpeg
├─ scripts
│  ├─ .DS_Store
│  ├─ .Rhistory
│  ├─ model_1
│  │  ├─ .DS_Store
│  │  ├─ model_1_testing.R
│  │  └─ model_1_trainign.R
│  ├─ model_2_6
│  │  ├─ .DS_Store
│  │  ├─ model_methylation_testing.R
│  │  ├─ model_methylation_training.R
│  │  ├─ model_mirna_testing.R
│  │  ├─ model_mirna_training.R
│  │  ├─ model_mrna_testing.R
│  │  ├─ model_mrna_training.R
│  │  ├─ model_mutation_testing.R
│  │  ├─ model_mutation_training.R
│  │  ├─ model_protein_testing.R
│  │  └─ model_protein_training.R
│  ├─ model_7
│  │  ├─ model_7_testing.R
│  │  └─ model_7_training.R
│  ├─ model_8
│  │  ├─ model_8_testing.R
│  │  └─ model_8_training.R
│  ├─ preprocessing.R
│  ├─ preprocessingFunctions.R
│  ├─ readingFiles.R
│  ├─ setup.R
│  └─ training.R
└─ sessionInfo.txt
```
Wokflow :
---------

data -> scripts -> results

Other files in the folder MLDM_assessment_2332635 :
---------------------------------------------------

1. environment.yml :
-----------------

This file contains all the dependencies and packages used in the entire project. The environment created by loading this file will act as a global environment containing all the dependecies specified.

2. The project is version controlled using git and the relevant details can be found in .git file. The .gitignore file is used to specift the files that need not be tracked by git and I have included history files and session data files in this.

3. The .snakemake file contains the logs and related data about the processing of snakemake pipelines.

4. The Snakefile contains the rules to run the snakemake pipeline.

Pipelining - Snakemake :
------------------------

The pipelining is done using snakemake pipelines (run the commands in activated conda environments):

1. The following command is used to run the pipeline and this creates the following folders and the files contained in it :
    
    rawData -> scripts -> results
    
  $ snakemake -cl


2. The clean rule is used to clean the directory of the files and folders generated by the snakemake pipeline :


  $ snakemake -cl clean

Git Version controlled :
------------------------


The following commands are used to version control using git :
(given git is already installed)

1. git init - to initialise a local git repository
2. git add . - to add all the files to be tracked by git
3. git commit -m "comments" - to commit the files and changes to staging area.
4. git git remote add origin https://github.com/your_username/your_repository.git
5. git push -u origin main
6. git clone /path/to/repository


Set up Snakemake pipeline :
---------------------------------

- Install conda, python, R and snakemake

- confirm the installation from the directory MLDM_assessment_2332635 :

    $which R 
    $which python
    $wnich conda
    
    The following commands will give the path where these are installed (if installations are run propoerly)

- conda might need to you to go to base directory (the case with me). Use the following command to open source in terminal :

      $source $HOME/miniconda3/bin/activate
      
      - or activate the conda base using :
      
          $conda activate base
      
- create environment.yml file and add the necessary dependency and sources

      $conda env create -f environment.yml

- activate the project conda environment

      $conda activate MLDM_assessment_2332635
    
- create a file Snakemake and write the rules inside. To run the snakeamke file first install and activate snakemake in the folder :

      $conda install -c conda-forge -c bioconda snakemake
      
- to deactivate the conda environment

      $conda activate MLDM_assessment_2332635

- once everything is run and the conda enironment is deactivated, you can activate the conda environment by simply running

      $conda activate MLDM_assessment_2332635
    
All the commands given above are to be run in the terminal in the project directory.

