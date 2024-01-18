#---
#title: "Model 7 training - omic gene and protein combined"
#author: 'StduentNo : 2332635'
#---

src1 <- paste(wd, "/scripts/training.R", sep = "")
source(src1)

#file_annotation_protein

## Read the protein annotation file
annotation_protein <- read.table(file_annotation_protein,header = TRUE)
head(annotation_protein)

distinct_values <- lapply(annotation_protein, unique)
distinct_values

# creating a new dataset and updating
protein_training <- training_data0$protein
protein_training$protein <- rownames(protein_training)

#View(protein_training)
# Read the gene data file
#omic_gene
gene_training <- training_data0$mrna
gene_training$gene <- rownames(gene_training)

#View(gene_training)

# Match proteins and genes based on protein annotation
matched_data_intermediate <- merge(protein_training, annotation_protein, by = "protein", all.x = TRUE, all.y = TRUE)

#View(matched_data_intermediate)
# Merge matched data with gene data based on gene
matched_data_protein_gene <- merge(matched_data_intermediate, gene_training, by = "gene")

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

#View(df_gene_x)
#View(df_protein_y)
# combining the two df's
combined_df_omics_training <- cbind(df_gene_x, df_protein_y)
combined_df_omics_training <- t(combined_df_omics_training)
#colnames(combined_df_omics_training)
#head(combined_df_omics_training)


class(combined_df_omics_training)
## remove features with no variance or all missing values

feature.var <- apply(combined_df_omics_training,1,var,na.rm=T)
combined_df_omics_training <- combined_df_omics_training[which(feature.var > 2e-16),]

class(combined_df_omics_training)
#View(combined_df_omics_training)
dim(combined_df_omics_training)
## identify top univariate predictors in each data type

## for each data type ...
cat(date(), "testing mrna and protein omics combined", "...\n")
## prepare to test, for each feature, feature ~ outcome
outcome <- training_clinical_numeric$pfi[match(colnames(combined_df_omics_training),rownames(training_clinical_numeric))]
design <- model.matrix(~outcome)
## fit linear model for each feature
fit <- lmFit(combined_df_omics_training, design)
## calculate p-values
fit <- eBayes(fit)
## identify the top 100 associations
idx <- order(fit$p.value[,"outcome"],decreasing=F)[1:25]
## return the names of the features with the top 25 associations
predictors_omics <- rownames(combined_df_omics_training)[idx]

nrow(combined_df_omics_training)
ncol(combined_df_omics_training)

## matching participants with data for clinical_numeric and omic_mirna

common_omics_ids_training <- rownames(training_clinical_numeric)
common_omics_ids_training <- intersect(common_omics_ids_training, colnames(combined_df_omics_training))

## construct a dataset with data on mirna (in transpose)
dataset_combined_omics_training <- t(combined_df_omics_training[predictors_omics,common_omics_ids_training])

# investigating the presence of missing data
d <- data.frame(dataset_combined_omics_training)
miss_var_summary(d)

# no missingness is observed in dataset_mirna_clinical as well as in the output, pfi
# as no missingness, imputation is excluded

## outcome variable
outcome_omics_training <- training_clinical_numeric[common_omics_ids_training,outcome.var]

outcome_omics_training <- as.factor(outcome_omics_training)

outcome_omics_training <- ifelse(outcome_omics_training==0,"No","Yes")

## standardize features (variance and scaling) 

dataset_combined_omics_training_var <- apply(dataset_combined_omics_training,2,var)
dataset_combined_omics_training <- dataset_combined_omics_training[,dataset_combined_omics_training_var > 2e-16]
dataset_combined_omics_training <- scale(dataset_combined_omics_training)

## specify controls and train the models
model_list_combined_omics <- train.model(dataset_combined_omics_training, outcome_omics_training)

########################################
# train_model_list
resamp <- resamples(model_list_combined_omics)

summary(resamp)

# plot the models

lattice::bwplot(resamp, metric = "ROC")
jpeg(file = "results/model_7/fig_roc_plot_combined_omics.jpeg")

dotplot(resamp)

# rf model
dataset_combined_omics_training <- cbind(dataset_combined_omics_training, outcome_omics_training)


train.predicts_combined_omics <- predict(
  model_list_combined_omics[["rf"]], 
  newdata = dataset_combined_omics_training)


cm <- confusionMatrix(train.predicts_combined_omics, as.factor(outcome_omics_training))

# Extract accuracy and 95% confidence interval
accuracy <- round(cm$overall["Accuracy"],digits = 3)
ci_lower <- round(cm$overall["AccuracyLower"],digits = 3)
ci_upper <- round(cm$overall["AccuracyUpper"],digits = 3)

# Add a new row to the data frame
model.accuracy.model_7 <- data.frame(Model = "rf", 
                                         Accuracy = accuracy, 
                                         CI_95_percent = paste(ci_lower, "-" , ci_upper)
)
