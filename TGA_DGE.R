# set working directory
setwd("~/Documents/smith_lab/TGA_3pod_analysis")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("tidyverse")
library(tidyverse)
BiocManager::install("GEOquery", force = TRUE)
library(GEOquery)
library(limma)
install.packages("edgeR")
library(edgeR)

# download geo dataset
geo_data <- getGEO("GSE215939", GSEMatrix = TRUE)
metadata <- pData(geo_data[[1]])

# Read raw data
files <- list.files(path = "./GSE215939_RAW", pattern = "\\.txt$", full.names = TRUE)
x <- read.maimages(files, source="agilent", green.only = TRUE, other.columns="gIsWellAboveBG")
colnames(x$E) <- gsub("^.*/(GSM[0-9]+)_.*$", "\\1", colnames(x))
colnames(x$Eb) <- gsub("^.*/(GSM[0-9]+)_.*$", "\\1", colnames(x))

# Create MA-plots for all the arrays
plotMA3by2(x)

# Create boxplots (pre-normalization)
library(ggplot2)
signal_data_frame <- stack(as.data.frame(log2(x$E)))
background_data_frame <- stack(as.data.frame(log2(x$Eb)))
ggplot(signal_data_frame) +
  geom_boxplot(aes(x = ind, y = values)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Signal")
ggplot(background_data_frame) +
  geom_boxplot(aes(x = ind, y = values))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Background")

y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")

Control <- y$genes$ControlType==1L
smallestGroupSize <- 5
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= smallestGroupSize
notProbe <- !grepl("^A_", y$genes$GeneName)
yfilt <- y[!Control & IsExpr &notProbe, ]

# Create boxplot for normalized data
signal_data_frame <- stack(as.data.frame(yfilt$E))
ggplot(signal_data_frame) +
  geom_boxplot(aes(x = ind, y = values)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Signal Normalized")

# Create dataframe for expression by patient, aggergate gene variants using mean
data_with_dupes = data.frame(GeneName = yfilt$genes$GeneName, Expression=yfilt$E)
data <- aggregate(. ~ GeneName, data = data_with_dupes, FUN = mean)
rownames(data) <- data$GeneName
data$GeneName <- NULL

colnames(data)<- gsub("^Expression.", "", colnames(data))
accession_to_subject <- setNames( metadata$title, metadata$geo_accession)
colnames(data) <- accession_to_subject[colnames(data)]

# Remove rows with negative values
#rows_with_negatives <- apply(data, 1, function(row) any(row < 0))
#data <- data[!rows_with_negatives, , drop = FALSE]

# Control_Male, Control_Female, TGALV_Male, TGALV_Female, TGARV_Male, TGARV_Female
metadata$source_name_ch1 <- gsub("-", "", metadata$source_name_ch1)
group <- factor(paste(metadata$source_name_ch1, metadata$`gender:ch1`, sep = "_"))
design <- model.matrix(~0 + group)
rownames(design) <- metadata$title
colnames(design) <- levels(group)

# Reorder design matrix to match column name in data
design <- design[match(colnames(data), rownames(design)), ]
all(rownames(design) == colnames(data))  # Should return TRUE

# Create DGE list
#group <- apply(design, 1, function(row) colnames(design)[which(row == 1)])

#dge_list <- DGEList(data, group = treatment_groups)

# Filter out lowly expressed genes
#keep <- filterByExpr(dge_list)
#dge_list <- dge_list[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
#dge_list <- normLibSizes(dge_list)
#dge_list$counts <- dge_list$counts*dge_list$samples$norm.factors

# sanity check with plotMDS
treatment_groups = unname(group)
group_mapping <- c("TGARV_Male" = 6, "TGARV_Female" = 5, "TGALV_Male" = 4, "TGALV_Female" = 3, "Control_Male" = 2, "Control_Female" = 1)
numeric_groups <- group_mapping[treatment_groups]
plotMDS(data, col = numeric_groups)

aw <- arrayWeights(data, design)

fit <- lmFit(data, design = design, weights = aw)

contrast_matrix <- makeContrasts(
  TGA_LV_Female_vs_Control_Female = TGALV_Female - Control_Female,
  TGA_LV_Male_vs_Control_Male = TGALV_Male - Control_Male,
  TGA_RV_Female_vs_Control_Female = TGARV_Female - Control_Female,
  TGA_RV_Male_vs_Control_Male = TGARV_Male - Control_Male,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = T)
results <- list()

for(i in colnames(contrast_matrix)){
  summary(fit2)
  results[[i]] <- topTable(fit2, coef = i, number = Inf) %>% rownames_to_column("Symbol")
}

# Iterate through 'results' and create a CSV for each entry
for(i in names(results)) {
  columns = results[[i]][, c("Symbol","logFC", "P.Value", "adj.P.Val")]
  write.csv(columns, file = paste0(i, ".csv"), row.names = FALSE)
}
