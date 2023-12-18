library(dplyr)
library(tidyverse)
library(GEOquery)

#read data
path = "GSE155237_gene_counts.csv"
data <- read.csv(path)
dim(data)


#get metadata
gse <- getGEO(GEO = 'GSE155237', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))


#metadata.subset <- select(metadata, c(1,13))
metadata.modified <- metadata %>%
  select(1,2,46)
#  rename(outcome = outcome:ch1) %>%
#  mutate(outcome = gsub("outcome: ", "", outcome)) %>%
 # mutate(outcome = gsub("Infected", "1", outcome)) %>%
#  mutate(outcome = gsub("Uninfected", "0", outcome))

head(data) 

#reshaping data to long format
#wide format, here all samples are each in a column
#long format, one column with genes, one with all samples, one with FPKM
#gather method convert from wide format to long format

data.long <- data %>%
  gather(key = 'samples', value='FPKM', -X) %>%
  mutate(samples = gsub("US.", "US-", samples))


#could be used for ML
#join dataframes = data.long + metadata.modified
data.long <- data.long %>%
  left_join(metadata.modified, by = c("samples" = "title"))



##############
####DESeq#####
##############

library(DESeq2)
library(tidyverse)
library(readr)

#Step 1: prepare counts data

#read data
path = "GSE155237_gene_counts.csv"
counts_data <- read.csv(path)

# make column 1 of counts data to become row names and delete column 1
row_names <- rownames(counts_data)
rownames(counts_data) <- counts_data$X
counts_data <- counts_data[, -1]

#read sample info
colData <- metadata.modified


colData <- colData %>%
  mutate(title = gsub("US-", "US.", title))

# rename colData from "outcome:ch1" to "Outcome"
colnames(colData)[colnames(colData) == "outcome:ch1"] <- "Outcome"


# replace column names in counts data with row names in colData
for (col_name in colnames(counts_data)) {
  matching_value <- colData$title[colData$title == col_name]
    if (length(matching_value) > 0) {
    replacement_value <- colData$geo_accession[colData$title == col_name]
    colnames(counts_data)[colnames(counts_data) == col_name] <- replacement_value
  }
}



#check row names in colData matches column names in counts_data 
all(colnames(counts_data) %in% rownames(colData))
#check row names in colData are in same order as column names in counts_data 
all(colnames(counts_data) == rownames(colData))
# Check for differences: row names in colData not in column names of counts_data
differences <- setdiff(rownames(colData), colnames(counts_data))


#Step 2: construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ Outcome)

keep <- rowSums(counts(dds) >= 35) >=6
dds <- dds[keep,]

#set factor level
dds$Outcome <- relevel(dds$Outcome, ref = "Uninfected")


#Step 3: Run DESeq
dds <- DESeq(dds)
resultsNames(dds)


############
res <- results(dds, alpha = 0.01)
res_infected_vs_uninfected <- results(dds, contrast = c("Outcome", "Infected", "Uninfected"), alpha = 0.01)
summary(res_infected_vs_uninfected)

res_infected <- results(dds, contrast = c("Outcome", "Infected", "Uninfected"))
DEGs_infected <- subset(res_infected, padj < 0.01 & log2FoldChange > 0.5)

res_uninfected <- results(dds, contrast = c("Outcome", "Uninfected", "Infected"))
DEGs_uninfected <- subset(res_uninfected, padj < 0.01 & log2FoldChange > 0.5)



##########

# Wald tests
res <- results(dds)

# Add Outcome information to the results DataFrame
outcome_column <- colData(dds)$Outcome
res$Outcome <- outcome_column

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

results <- results(dds, contrast = c("Outcome", "Infected", "Uninfected"))
results


# MA plot
plotMA(res)




