#---
#title: "Model 7 testing - omic gene and protein combined"
#author: 'StduentNo : 2332635'
#---


#file_annotation_protein

head(annotation_protein)


# creating a new dataset and updating
protein_testing <- testing_data0$protein
protein_testing$protein <- rownames(protein_testing)

#View(protein_testing)
# Read the gene data file
#omic_gene
gene_testing <- testing_data0$mrna
gene_testing$gene <- rownames(gene_testing)

#View(gene_testing)

# Match proteins and genes based on protein annotation
matched_data_intermediate <- merge(protein_testing, annotation_protein, by = "protein", all.x = TRUE, all.y = TRUE)

#View(matched_data_intermediate)
# Merge matched data with gene data based on gene
matched_data_protein_gene <- merge(matched_data_intermediate, gene_testing, by = "gene")

# protein-gene pair
protein_gene_pair <- matched_data_protein_gene[,c("gene","protein")]

# rows - patient ids
# columns - each gene-protein combo.

# new df for sample.x - protein
df_gene_x <- matched_data_protein_gene[,c("gene","protein")]

# Loop over column names and add all .x to new data frame
for (col in names(matched_data_protein_gene)) {
  if (grepl("\\.x$", col)) {
    sample <- sub("\\.x$", "", col)
    df_gene_x[[sample]] <- matched_data_protein_gene[, col]
  } 
}

colnames(df_gene_x)
head(df_gene_x)


# replacing all the symbols with "_"
df_gene_x$protein <- gsub("[^[:alnum:]_]+", "_", df_gene_x$protein)
df_gene_x$gene <- gsub("[^[:alnum:]_]+", "_", df_gene_x$gene)

df_gene_x$protein_gene_pair <- paste(df_gene_x$protein, df_gene_x$gene, sep = "_._")


df_gene_x <- df_gene_x %>% select(-c("gene","protein"))
df_gene_x <- data.frame(protein_gene_pair = df_gene_x$protein_gene_pair, df_gene_x[, -which(names(df_gene_x) == "protein_gene_pair")])

rownames(df_gene_x) <- df_gene_x$protein_gene_pair
df_gene_x <- df_gene_x %>% select(-c("protein_gene_pair"))


df_gene_x <- data.frame(t(df_gene_x))

# specifying if each column represents gene or protein expression 
colnames(df_gene_x) <- paste0(colnames(df_gene_x), "_._gene")

df_gene_x

# new df for sample.y - gene

df_protein_y <- matched_data_protein_gene[,c("gene","protein")]

# Loop over column names and add all .x to new data frame
for (col in names(matched_data_protein_gene)) {
  if (grepl("\\.y$", col)) {
    sample <- sub("\\.y$", "", col)
    df_protein_y[[sample]] <- matched_data_protein_gene[, col]
  } 
}

colnames(df_protein_y)
head(df_protein_y)


# replacing all the symbols with "_"
df_protein_y$protein <- gsub("[^[:alnum:]_]+", "_", df_protein_y$protein)
df_protein_y$gene <- gsub("[^[:alnum:]_]+", "_", df_protein_y$gene)

df_protein_y$protein_gene_pair <- paste(df_protein_y$protein, df_protein_y$gene, sep = "_._")


df_protein_y <- df_protein_y %>% select(-c("gene","protein"))
df_protein_y <- data.frame(protein_gene_pair = df_protein_y$protein_gene_pair, df_protein_y[, -which(names(df_protein_y) == "protein_gene_pair")])

rownames(df_protein_y) <- df_protein_y$protein_gene_pair
df_protein_y <- df_protein_y %>% select(-c("protein_gene_pair"))


df_protein_y <- data.frame(t(df_protein_y))

# specifying if each column represents gene or protein expression 
colnames(df_protein_y) <- paste0(colnames(df_protein_y), "_._protein")

df_protein_y

View(df_gene_x)
View(df_protein_y)
# combining the two df's
combined_df_omics_testing <- cbind(df_gene_x, df_protein_y)
combined_df_omics_testing <- t(combined_df_omics_testing)
#colnames(combined_df_omics_testing)
#head(combined_df_omics_testing)


class(combined_df_omics_testing)
## remove features with no variance or all missing values

feature.var <- apply(combined_df_omics_testing,1,var,na.rm=T)
combined_df_omics_testing <- combined_df_omics_testing[which(feature.var > 2e-16),]

class(combined_df_omics_testing)
#View(combined_df_omics_testing)
dim(combined_df_omics_testing)
## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna and protein omics combined", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- testing_clinical_numeric$pfi[match(colnames(combined_df_omics_testing),rownames(testing_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(combined_df_omics_testing, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omics <- rownames(combined_df_omics_testing)[idx]

nrow(combined_df_omics_testing)
ncol(combined_df_omics_testing)

## matching participants with data for clinical_numeric and omic_mirna

common_omics_ids_testing <- rownames(testing_clinical_numeric)
common_omics_ids_testing <- intersect(common_omics_ids_testing, colnames(combined_df_omics_testing))

## construct a dataset with data on mirna (in transpose)
dataset_combined_omics_testing <- t(combined_df_omics_testing[predictors_omics,common_omics_ids_testing])

# investigating the presence of missing data
d <- data.frame(dataset_combined_omics_testing)
miss_var_summary(d)

# no missingness is observed in dataset_mirna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_omics_testing <- testing_clinical_numeric[common_omics_ids_testing,outcome.var]

outcome_omics_testing <- as.factor(outcome_omics_testing)

outcome_omics_testing <- ifelse(outcome_omics_testing==0,"No","Yes")

## standardize features (variance and scaling) 

dataset_combined_omics_testing_var <- apply(dataset_combined_omics_testing,2,var)
dataset_combined_omics_testing <- dataset_combined_omics_testing[,dataset_combined_omics_testing_var > 2e-16]
dataset_combined_omics_testing <- scale(dataset_combined_omics_testing)

# predic

dataset_combined_omics_testing <- cbind(dataset_combined_omics_testing, outcome_omics_testing)


train.predicts_combined_omics <- predict(
  model_list_combined_omics[["rf"]], 
  newdata = dataset_combined_omics_testing)

confusionMatrix(train.predicts_combined_omics, as.factor(outcome_omics_testing))


