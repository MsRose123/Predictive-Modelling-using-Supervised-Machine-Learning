rule all:
    "The default rule"
    input:
        "results/model_1/fig_roc_plot_clinical.jpeg",
        "results/model_2_6/fig_roc_plot_mirna.jpeg",
        "results/model_2_6/fig_roc_plot_mrna.jpeg",
        "results/model_2_6/fig_roc_plot_mutation.jpeg",
        "results/model_2_6/fig_roc_plot_methylation.jpeg",
        "results/model_2_6/fig_roc_plot_protein.jpeg",
        "results/model_2_6/fig_roc_plot_combined_omics.jpeg",
        "results/model_2_6/fig_roc_plot_combined_clinical_omics.jpeg",
        

rule setup:
    "setting up, reading directory, loading data and preprocessing functions"
    shell: """
    Rscript scripts/setup.R,
    Rscript scripts/preprocessingFunctions.R,
    Rscript scripts/preprocessing.R,
    Rscript scripts/readingFiles.R,
    Rscript scripts/training.R
    """

rule model_1:
    "model 1 - Clinical features"
    shell: """
    mkdir -p cleanData
    Rscript scripts/model_1/model_1_trainign.R,
    Rscript scripts/model_1/model_1_testing.R
    """

rule model_2_6:
    "model 2 to 6 - multi-omics features"
    shell: """
    Rscript scripts/model_2_6/model_mirna_training.R,
    Rscript scripts/model_2_6/model_mrna_training.R,
    Rscript scripts/model_2_6/model_protein_training.R,
    Rscript scripts/model_2_6/model_mutation_training.R,
    Rscript scripts/model_2_6/model_methylation_training.R,
    Rscript scripts/model_2_6/model_mirna_testing.R,
    Rscript scripts/model_2_6/model_mrna_testing.R,
    Rscript scripts/model_2_6/model_protein_testing.R,
    Rscript scripts/model_2_6/model_mutation_testing.R,
    Rscript scripts/model_2_6/model_methylation_testing.R
    """

rule model_7:
    "model 1 - protein and gene features"
    shell: """
    mkdir -p cleanData
    Rscript scripts/model_7/model_7_training.R,
    Rscript scripts/model_7/model_8_testing.R
    """

rule model_8:
    "model 1 - protein and gene combined with clinical features"
    shell: """
    mkdir -p cleanData
    Rscript scripts/model_8/model_8_training.R,
    Rscript scripts/model_8/model_8_testing.R
    """


rule clean:
    "Clean up"
    shell: """
    if [ -d cleanData ]; then
      rm -r cleanData
    else
      echo directory cleanData does not exist
    fi
    if [ -d intermediates ]; then
      rm -r intermediates
    else
      echo directory intermediates does not exist
    fi
    if [ -d results ]; then
      rm -r results
    else
      echo directory results does not exist
    fi
    """
