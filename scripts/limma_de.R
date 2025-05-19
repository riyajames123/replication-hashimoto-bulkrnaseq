# Differential Expression Analysis using limma + voom
# Author: Riya James (with ChatGPT)
# Description: This script performs differential expression analysis on RNA-seq count data using limma + voom

# ===================
# 1. Load Libraries
# ===================
library(edgeR)      # for DGEList, normalization
library(limma)      # for voom and differential testing

# ============================
# 2. Load and Prepare Data
# ============================
# Replace this with your actual file path
data_file <- "/scratch/james.ri/trans_proj/featurecounts_new_output/gene_counts.tsv"

# Read raw counts (skip comment lines)
raw_counts <- read.table(data_file, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)

# Remove unnecessary columns (Chr, Start, End, etc.)
# Keeping only Geneid and count columns
count_data <- raw_counts[, c(1, 7:ncol(raw_counts))]

# Set Gene IDs as rownames
rownames(count_data) <- count_data$Geneid
count_data <- count_data[, -1]  # Drop Geneid column now that it's rownames

# clean up the column names
colnames(count_data) <- gsub(".*/|\\.sorted\\.bam", "", colnames(count_data))


# ============================
# 3. Create Sample Metadata
# ============================
sample_info <- data.frame(
    sample = c("HRR568836", "HRR568837", "HRR568838", "HRR568839", "HRR568840",
               "HRR568844", "HRR568849", "HRR568850", "HRR568852", "HRR568857"),
    condition = factor(c("control", "control", "control", "control", "control",
                         "disease", "disease", "disease", "disease", "disease")),
    sex = factor(c("female", "female", "male", "male", "female",
                   "female", "female", "female", "female", "female"))
)
rownames(sample_info) <- sample_info$sample


# ============================
# 4. Create DGEList Object
# ============================
dge <- DGEList(counts = count_data, group = sample_info$condition)
dge <- calcNormFactors(dge)  # Normalize using TMM

# ============================
# 5. voom Transformation
# ============================
design <- model.matrix(~ condition, data = sample_info)
voom_data <- voom(dge, design, plot = TRUE)  # Generates mean-variance plot

# ============================
# 6. Fit Linear Model
# ============================
fit <- lmFit(voom_data, design)
fit <- eBayes(fit)

# ============================
# 7. Extract DE Results
# ============================
# coef=2 corresponds to conditiondisease if "control" is baseline
results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Save to file
write.table(results, "DE_results_limma_voom.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# Optional: View top hits
head(results)

# ============================
# 8. Volcano Plot (optional)
# ============================
library(ggplot2)
results$gene <- rownames(results)
results$significant <- results$adj.P.Val < 0.05

ggplot(results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("black", "red")) +
    labs(title = "Differential Gene Expression- Disease vs Control",
         x = "log2 Fold Change", y = "-log10 p-value")
