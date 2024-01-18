#!/bin/bash

#source $HOME/miniconda3/bin/activate

# Enter the base conda environment and make sure it is up to date
conda activate base
conda update -n base -c defaults conda

# Create the conda environment
conda env create -f environment.yml

# To activate the conda environment
conda activate MLDM_assessment_2332635

#clean the pipeline
snakemake -c1 clean

# run pipeline
snakemake -c1

