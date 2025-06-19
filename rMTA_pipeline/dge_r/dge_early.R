# -----------------------------
# Load necessary libraries
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GEOquery", "edgeR", "limma", "ggplot2", "reshape2"))


library(GEOquery)
library(edgeR)
library(limma)
library(ggplot2)
library(reshape2)

# Download supplementary files for GSE203206 [https://molecularbrain.biomedcentral.com/articles/10.1186/s13041-022-00963-2#MOESM2] if the directory doesn't exist.
if (!dir.exists("GSE203206")) {
    getGEOSuppFiles("GSE203206", makeDirectory = TRUE, baseDir = ".")
}


# -----------------------------
# Part I: Differential Expression Analysis Using limma-voom
# -----------------------------

# Step 1: Read in the gzipped counts and metadata files from
counts <- read.table("./GSE203206/GSE203206_Subramaniam.ADRC_brain.counts.tsv.gz",
                     header = TRUE, sep = "\t", row.names = 1)

# Read in the metadata file
metadata <- read.table("./GSE203206/GSE203206_Subramaniam.ADRC_brain.metadata.tsv.gz",
                       header = TRUE, sep = "\t", row.names = 1)


# Preview the first few rows of counts and metadata
cat("Counts matrix (first few rows):\n")
print(head(counts))
cat("\nMetadata (first few rows):\n")
print(head(metadata))

# Step 2: Ensure that the sample names match between counts and metadata.
# Here, we reorder metadata to match the counts columns.
metadata <- metadata[colnames(counts), ]

metadata$Sample <- substr(metadata$Sample, 1, 3)

# if Sample starts with "LOL", remove that line
metadata <- metadata[!grepl("^LOL", metadata$Sample), ]

# remove NAs from metadata
metadata <- metadata[!is.na(metadata$Sample), ]

# remove from counts the columns that are not in metadata
counts <- counts[, colnames(counts) %in% rownames(metadata)]




# Step 3: Create a DGEList object from the counts data
dge <- DGEList(counts = counts)

# Step 4: Normalize the data using TMM normalization
dge <- calcNormFactors(dge)


# Just keep the first 3 letter in the Sample column of metadata




# Step 5: Prepare the design matrix using a grouping variable.
# Here we assume 'PATHDX' (pathological diagnosis) defines the groups.
metadata$Sample <- factor(metadata$Sample)
design <- model.matrix(~ Sample, data = metadata)
cat("\nDesign matrix (first few rows):\n")
print(head(design))

# Step 6: Apply the voom transformation to model the mean-variance relationship
v <- voom(dge, design, plot = TRUE)

# Step 7: Fit the linear model and apply empirical Bayes moderation
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Step 8: Extract differential expression results.
# coef = 2 corresponds to the effect of PATHDX (assuming two groups, e.g., Alzheimer’s vs. Control)
results_limma <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

# Add the gene identifier column
results_limma$gene_id <- rownames(results_limma)
final_output <- results_limma[, c("gene_id", "logFC", "adj.P.Val")]

# Save the DGE results to a CSV file
write.csv(final_output, file = "DGE_results_limma_voom_EOL.csv", row.names = FALSE)
cat("\nDifferential Expression Analysis Complete. Results saved as 'DGE_results_limma_voom.csv'.\n")

# -----------------------------
# Part II: Metadata Analysis
# -----------------------------

cat("\n------ Metadata Analysis ------\n")

# Display summary statistics for all metadata columns
cat("\nSummary of Metadata:\n")
print(summary(metadata))

# Frequency tables for categorical variables
cat("\nFrequency Table for PATHDX:\n")
print(table(metadata$PATHDX))

cat("\nFrequency Table for SEX:\n")
print(table(metadata$SEX))

# Descriptive statistics for numeric variables (e.g., AGE and EDUCATION)
cat("\nDescriptive Statistics for AGE:\n")
cat("Mean Age:", mean(metadata$AGE, na.rm = TRUE), "\n")
cat("Median Age:", median(metadata$AGE, na.rm = TRUE), "\n")
cat("Age Range:", paste(range(metadata$AGE, na.rm = TRUE), collapse = " - "), "\n")

cat("\nDescriptive Statistics for EDUCATION:\n")
cat("Mean Education:", mean(metadata$EDUCATION, na.rm = TRUE), "\n")
cat("Median Education:", median(metadata$EDUCATION, na.rm = TRUE), "\n")
cat("Education Range:", paste(range(metadata$EDUCATION, na.rm = TRUE), collapse = " - "), "\n")


cat("\nMetadata analysis complete.\n")
