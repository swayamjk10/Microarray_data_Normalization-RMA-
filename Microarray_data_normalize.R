install.packages("BiocManager") # Install BiocManager if not already installed
BiocManager::install("affy")    # Use BiocManager to install the affy package


library(affy)
library(GEOquery)
library(tidyverse)

# get supplementary files
getGEOSuppFiles("GSE148537")

# untar files
untar("GSE148537/GSE148537_RAW.tar", exdir = 'microarraydata/')

# reading in .cel files
raw_data <- ReadAffy(celfile.path = "microarraydata/")

# performing RMA normalization
rma_normalized_data <- rma(raw_data)

# get expression estimates
normalized_expr_data <- as.data.frame(exprs(rma_normalized_data))


# map probe IDs to gene symbols
gse <- getGEO("GSE148537", GSEMatrix = TRUE)


# fetch feature data to get ID - gene symbol mapping
feature_data <- gse$GSE148537_series_matrix.txt.gz@featureData@data

feature_data <- feature_data[,c(1,11)]

# joining feature data with normalised expression data
normalized_expr_data <- normalized_expr_data %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature_data, by = 'ID')


# Extract all unique genes from the column
all_genes <- unique(normalized_expr_data$ID) 


