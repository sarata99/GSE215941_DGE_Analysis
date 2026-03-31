# Set working directory and load packages
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


#------------------------ mRNA analysis ------------------------#
# Download geo dataset
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

# Control_Male, Control_Female, TGALV_Male, TGALV_Female, TGARV_Male, TGARV_Female
metadata$source_name_ch1 <- gsub("-", "", metadata$source_name_ch1)
group <- factor(paste(metadata$source_name_ch1, metadata$`gender:ch1`, sep = "_"))
design <- model.matrix(~0 + group)
rownames(design) <- metadata$title
colnames(design) <- levels(group)

# Reorder design matrix to match column name in data
design <- design[match(colnames(data), rownames(design)), ]
all(rownames(design) == colnames(data))  # Should return TRUE

# Sanity check with plotMDS
treatment_groups = unname(group)
group_mapping <- c("TGARV_Male" = 6, "TGARV_Female" = 5, "TGALV_Male" = 4, "TGALV_Female" = 3, "Control_Male" = 2, "Control_Female" = 1)
numeric_groups <- group_mapping[treatment_groups]
plotMDS(data, col = numeric_groups)

# Fix heteroscedasticity
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



#------------------------ miRNA analysis ------------------------#
geo_data <- getGEO("GSE215940", GSEMatrix = TRUE)
metadata <- pData(geo_data[[1]])

# Read raw data
files <- list.files(path = "./GSE215940_RAW", pattern = "\\.txt$", full.names = TRUE)
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

# Control_Male, Control_Female, TGALV_Male, TGALV_Female, TGARV_Male, TGARV_Female
metadata$source_name_ch1 <- gsub("-", "", metadata$source_name_ch1)
group <- factor(paste(metadata$source_name_ch1, metadata$`gender:ch1`, sep = "_"))
design <- model.matrix(~0 + group)
rownames(design) <- metadata$title
colnames(design) <- levels(group)

# Reorder design matrix to match column name in data
design <- design[match(colnames(data), rownames(design)), ]
all(rownames(design) == colnames(data))  # Should return TRUE

# Sanity check with plotMDS
treatment_groups = unname(group)
group_mapping <- c("TGARV_Male" = 6, "TGARV_Female" = 5, "TGALV_Male" = 4, "TGALV_Female" = 3, "Control_Male" = 2, "Control_Female" = 1)
numeric_groups <- group_mapping[treatment_groups]
plotMDS(data, col = numeric_groups)

# Fix heteroscedasticity
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
  write.csv(columns, file = paste0(i, "_miRNA.csv"), row.names = FALSE)
}


#------------------------ mi-RNA-mRNA interaction ------------------------#
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
library(httr)
library(tidyverse)

#--- Official mirDIP API helper functions ---#

url <- "http://ophid.utoronto.ca/mirDIP"

mapScore <- list("0", "1", "2", "3")
names(mapScore) <- c("Very High", "High", "Medium", "Low")

unidirectionalSearchOnMicroRNAs <- function(microRNAs, minimumScore) {
  parameters <- list(
    genesymbol = "",
    microrna = microRNAs,
    scoreClass = unlist(mapScore[minimumScore]) # Unlist ensures pure string for POST
  )
  # Send http POST to the specific Unidirectional endpoint
  res <- POST(paste(url, "/Http_U", sep = ""), body = parameters, encode = "form")
  return(res)
}

makeMap <- function(res) {
  ENTRY_DEL = "\001"
  KEY_DEL = "\002"
  response = content(res, "text", encoding = "UTF-8")
  arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))
  
  list_map <- list()
  vec_map_names <- character()
  
  for (str in arr) {
    arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE))
    if (length(arrKeyValue) > 1) {
      list_map <- c(list_map, list(arrKeyValue[2]))
      vec_map_names <- c(vec_map_names, arrKeyValue[1])
    }
  }
  names(list_map) <- vec_map_names
  return(list_map)
}

# For each comparison (TGA_LV_Female_vs_Control_Female, TGA_RV_Male_vs_Control_Male, etc)
for(i in colnames(contrast_matrix)) {
  # Grab mRNA and miRNA DEG files
  mRNA_file <- paste0(i, ".csv")
  miRNA_file <- paste0(i, "_miRNA.csv")
  
  if(file.exists(mRNA_file) & file.exists(miRNA_file)){
    # Can use adjusted p-value if want to be more strict
    mRNA_res <- read.csv(mRNA_file) %>% filter(P.Value < 0.05)
    miRNA_res <- read.csv(miRNA_file) %>% filter(P.Value < 0.05)
  }
  
  if(nrow(miRNA_res) == 0) {
    message(paste("No significant miRNAs found for", i, ". Skip miRNA query"))
    next
  }
  
  message(paste("Calling mirDIP API for", i))
  
  # API needs a string of miRNA's separated by commas
  mirna_string <- paste(miRNA_res$Symbol, collapse = ",")
  
  # Call API
  response <- unidirectionalSearchOnMicroRNAs(microRNAs = mirna_string, minimumScore = "Very High")
  
  if(status_code(response) == 200) {
    
    # Use the official makeMap function to unpack the envelope
    list_map <- makeMap(response)
    
    # Extract the raw table data from the 'results' key
    raw_table_data <- list_map[["results"]]
    results_size <- as.numeric(list_map[["results_size"]])
    
    if(is.null(raw_table_data) || is.na(results_size) || results_size == 0) {
      message(paste("API returned no targets for", i, "Skip miRNA query"))
      next
    }
    
    # Parse the formatted tab-delimited spreadsheet string into a dataframe
    target_df <- read_tsv(raw_table_data, show_col_types = FALSE) %>%
      rename(
        GeneSymbol = `Gene Symbol`,
        IntegratedScore = `Integrated Score`,
        NumberOfSources = `Number of Sources`
      ) %>%
      select(MicroRNA, GeneSymbol, IntegratedScore, NumberOfSources) %>%
      distinct()
    
  } else {
    message(paste("API Error: Server returned status", status_code(response), "\n"))
    next 
  }
  
  integrated_df <- target_df %>%
    inner_join(miRNA_res, by = c("MicroRNA" = "Symbol")) %>%
    rename(miRNA_logFC = logFC, miRNA_P = P.Value) %>%
    select(MicroRNA, GeneSymbol, IntegratedScore, NumberOfSources, miRNA_logFC, miRNA_P)
  
  integrated_df <- integrated_df %>%
    inner_join(mRNA_res, by = c("GeneSymbol" = "Symbol")) %>%
    rename(mRNA_logFC = logFC, mRNA_P = P.Value) %>%
    select(MicroRNA, GeneSymbol, IntegratedScore, NumberOfSources, miRNA_logFC, mRNA_logFC, miRNA_P, mRNA_P)
  
  # Filter for Inverse Expression 
  final_network <- integrated_df %>%
    filter((miRNA_logFC > 0 & mRNA_logFC < 0) | (miRNA_logFC < 0 & mRNA_logFC > 0)) %>%
    arrange(desc(IntegratedScore)) 
  
  # Save the final mapped network
  out_file <- paste0(i, "_mirDIP_API_Network.csv")
  write.csv(final_network, out_file, row.names = FALSE)
  
  message(paste("Saved", nrow(final_network), "validated inverse interactions to", out_file, "\n"))
  
}
  
