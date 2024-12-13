comparing transcriptome results mapped from different reference genome
================

## Intro

Import read count data from two files, with one (9192_counts.txt) mapped
into whole genome (hg19), and the other file
(9192_mini_gene_id_name_counts.txt) mapped to total transcriptome
reference .

## Read in data

``` r
# Load necessary library
library(tidyr)

# Read the file 
file_path <- "raw_data/9192_mini_gene_id_name_counts.txt"
data <- read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE)

# Split the second column into gene_id and gene_name
data <- separate(data, V2, into = c("gene_id", "gene_name"), sep = "\\|")

# Set appropriate column names
colnames(data) <- c("count", "gene_id", "gene_name")

# Preview the data
head(data)
```

    ##   count            gene_id gene_name
    ## 1  6394 ENSG00000000003.16    TSPAN6
    ## 2     2  ENSG00000000005.6      TNMD
    ## 3  5550 ENSG00000000419.14      DPM1
    ## 4  1710 ENSG00000000457.14     SCYL3
    ## 5  3024 ENSG00000000460.17     FIRRM
    ## 6   472 ENSG00000000971.17       CFH

``` r
# Check if there are duplicate gene names
any(duplicated(data$gene_name))
```

    ## [1] TRUE

``` r
# # If TRUE, print out the duplicates
# duplicates <- data[duplicated(data$gene_name) | duplicated(data$gene_name, fromLast = TRUE), ]
# 
# # View the duplicated gene names
# duplicates

any(duplicated(data$gene_id))
```

    ## [1] FALSE

``` r
##proceed analysis with ensembl ID
# Remove the version suffix from Ensembl IDs
data$gene_id_clean <- sub("\\..*", "", data$gene_id)

# # View the updated data
# head(data)
```

``` r
# Load necessary library
library(tidyr)

# Read the file 
file_path2 <- "raw_data/9192_counts.txt"
data_hg <- read.table(file_path2, header = TRUE, sep = "", stringsAsFactors = FALSE)
```

## Compare and visulize the read mapping results in both file

``` r
#colnames(data_hg)

# Load necessary libraries
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)

# Select only relevant columns before merging
data_selected <- data %>%
  select(gene_id_clean, count)  # Select 'gene_id_clean' and 'count' from 'data'

data_hg_selected <- data_hg %>%
  select(Geneid, `X9192_paired_sorted.bam`)  # Select 'Geneid' and the corresponding column from 'data_hg'

# Merge the two datasets by the cleaned gene ID
merged_data <- merge(data_selected, data_hg_selected, by.x = "gene_id_clean", by.y = "Geneid")

# Rename the columns for easier reference
colnames(merged_data) <- c("gene_id", "count_data", "count_data_hg")


merged_data_log <- merged_data

# Apply log(read + 1) transformation to both read counts to handle zeros
merged_data_log$count_data <- log10(merged_data$count_data + 1)
merged_data_log$count_data_hg <- log10(merged_data$count_data_hg + 1)




# Create a dot plot using ggplot2 with transformed counts
ggplot(merged_data_log, aes(x = count_data, y = count_data_hg)) +
  geom_point() +  # Scatter plot
  labs(
    title = "Comparison of Read Counts Between Two Datasets (Log(Count + 1) Scale)",
    x = "Log(Count + 1) from First Table (data)",
    y = "Log(Count + 1) from Second Table (data_hg)"
  ) +
  theme_minimal()  # Clean theme for visualization
```

![](compare_transcriptome_mapped_files/figure-gfm/compare%20results-1.png)<!-- -->

``` r
# Calculate Pearson correlation for log-transformed counts
cor_log <- cor(merged_data_log$count_data, merged_data_log$count_data_hg, method = "pearson")

# Print the correlation value for log-transformed counts
cat("Pearson Correlation (log-transformed counts):", cor_log, "\n")
```

    ## Pearson Correlation (log-transformed counts): 0.8289685

## Conclusion

Acceptable correlation was observed for the table mapped from different
references. Proceed with mapping to total transcriptome reference, as it
is more time efficient when mapping.
