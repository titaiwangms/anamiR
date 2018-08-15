# anamiR-2.0
An integrated analysis R package of miRNA and mRNA profiiling

## anamiR

This document guides the user through all available functions of the anamiR package. anamiR aims to find potential miRNA-target gene interactions from both miRNA and mRNA expression data.

Traditional miRNA analysis methods use online databases to predict miRNA-target gene interactions. However, the inconsistent results make the interactions less reliable. To address this issue, anamiR integrates the whole expression analysis with expression data into a workflow including normalization, differential expression, correlation, and finally database intersection, to find more reliable interactions.
Moreover, users can identify interactions from pathways or gene sets of interest.

## Data Source

As shown in the workflow, not only samples of *paired* miRNA and mRNA expression data, but also phenotypic information on miRNA and mRNA sequences are required for the analysis. Since anamiR reads data in expression matrices, data sources are platform- and technology-independent. Particularly, expression data from microarrays or next generation sequencing are all acceptable for anamiR. However, this also means the raw data have to be formatted into expression matrices before using anamiR.

### mRNA expression

Columns for samples. Rows for genes

```raw
GENE  SmapleA   SamplB ...
A     0.1       0.2
B    -0.5      -0.3
C     0.4       0.1
```

### miRNA expression

Columns for samples. Rows for miRNAs

```raw
miRNA  SampleA  SampleB ...
A         0.1     0.5
B        -0.5     2.1
C         0.4     0.3
```

### phenotype data

Cloumns for samples. Rows for feature name, including two groups, multiple groups, continuous data.

```raw
Feature  groupA  groupA  groupA ...
SampleA  123.33     A       A
SampleB  120.34     B       C
SampleC  121.22     A       B
```


## Installation

anamiR is on Bioconductor and can be installed following standard installation procedure.

```R
source("http://www.bioconductor.org/biocLite.R")
biocLite("anamiR")
```

To use,

```R
library(anamiR)
```

## General Workflow

Basically there are six steps, corresponding to six R functions, to complete the whole analysis:

1. Normalize expression data
2. Find differentially expressed miRNAs and genes
3. Convert miRNA annotation to the latest version
4. Calculate the correlation coefficient between each miRNA and gene
5. Predict and validate the intersection of miRNA-target gene interaction databases
6. Functional analysis of genes of interest

![alt text](https://github.com/AllenTiTaiWang/anamiR-2.0/blob/master/vignettes/pics/Generalworkflow.png)

1. Import data
```R
library(anamiR)
data(mrna)
data(mirna)
data(pheno.mirna)
data(pheno.mrna)
```

2. SummarizedExperimaent object

Before entering the main workflow, we should put our data and phenotype 
information into “SummarizedExperiment” format first.

```R
mirna_se <- SummarizedExperiment(assays = SimpleList(counts=mirna), colData = pheno.mirna)
mrna_se <- SummarizedExperiment(assays = SimpleList(counts=mrna), colData = pheno.mrna)
```

3. Differential Expression analysis

Second, we will find differentially expressed genes and miRNAs. There are three statistical methods in this function. Here, we use “limma” for demonstration.

```R
mirna_d <- differExp_discrete(se = mirna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
mrna_d <- differExp_discrete(se = mrna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
```

This function will delete genes and miRNAs (rows) that are not differentially expressed, and add another three columns to represent fold-change (log2), p-value, adjusted p-value.

4. miRNA ID conversion

Before using the collected databases to find the intersection of potential miRNA-target gene interactions, we have to make sure all miRNA data are in the latest annotation version (miRBase 21). If not, we can use this step to do it.

```R
mirna_21 <- miR_converter(data = mirna_d, remove_old = TRUE, original_version = 17, latest_version = 21)
```

5. Correlation abalysis

Fourth, to find potential miRNA-target gene interactions, we have to combine the information from two differentially expressed datasets, which we obtain from “differExp_discrete”.

```R
cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_21, cut.off = -0.5, method = "pearson")
```

In the displayed list, each row is a potential interaction, and only the rows whose correlation coefficient is < cut.off would be kept in the list.

Note that in our assumption, miRNAs negatively regulate the expression of their target genes; that is, cut.off should typically be a negative decimal.

6. Database Intersection

After correlation analysis, we have some potential interactions, and we use the “database_support” function to get information on whether there are databases that predict or validate these interactions.

```R
sup <- database_support(cor_data = cor, Sum.cutoff = 2, org = "hsa")
```

In the output, the Sum column tells us the total hits in the 8 prediction databases and the Validate column tells us if each interaction has been validated.

7. Heatmap visualization

```R
heat_vis(cor, mrna_d, mirna_21)
```

8. Functional analysis

Finally, after finding reliable miRNA-target gene interactions, we are also interested in the pathways that may be enriched by the expression of these genes.

```R
enr <- enrichment(data_support = sup, per_time = 5000)
```

Note that for the parameter “per_time”, we chose 500 times for demonstration purposes. The default value is 5000 times.

The output from this data shows not only the P-value generated by the hypergeometric test, but also the empirical P-value, which is the value of the average permutation test in each pathway.

## GSEA Workflow

Basically there are only two steps with two R functions, to complete the whole analysis:

1. Find related miRNAs and genes in the possible enriched pathways.
2. Find potential interactions from the above result.

![alt text](https://github.com/AllenTiTaiWang/anamiR-2.0/blob/master/vignettes/pics/GSEAworkflow.png)

1. Load data

```R
aa <- system.file("extdata", "GSE19536_mrna.csv", package = "anamiR")
mrna <- data.table::fread(aa, fill = T, header = T)
bb <- system.file("extdata", "GSE19536_mirna.csv", package = "anamiR")
mirna <- data.table::fread(bb, fill = T, header = T)
cc <- system.file("extdata", "pheno_data.csv", package = "anamiR")
pheno.data <- data.table::fread(cc, fill = T, header = T)
```

2. transform data format

```R
mirna_name <- mirna[["miRNA"]]
mrna_name <- mrna[["Gene"]]
mirna <- mirna[, -1]
mrna <- mrna[, -1]
mirna <- data.matrix(mirna)
mrna <- data.matrix(mrna)
row.names(mirna) <- mirna_name
row.names(mrna) <- mrna_name
pheno_name <- pheno.data[["Sample"]]
pheno.data <- pheno.data[, -1]
pheno.data <- as.matrix(pheno.data)
row.names(pheno.data) <- pheno_name
```

#### Extra steps

```R
mirna <- mirna[, 50:70]
mrna <- mrna[, 50:70]
pheno.data <- as.matrix(pheno.data[50:70, ])
colnames(pheno.data) <- "ER"
```

3. SummarizedExperiment object

```R
mirna_21 <- miR_converter(mirna, original_version = 17)
mirna_se <- SummarizedExperiment(assays = SimpleList(counts=mirna_21), colData = pheno.data)
mrna_se <- SummarizedExperiment(assays = SimpleList(counts=mrna), colData = pheno.data)
```

4. GSEA_ana

In the first step, we use the “GSEA_ana” function to find the pathways which are the most likely to be enriched in the given expression data.

```R
table <- GSEA_ana(mrna_se = mrna_se, mirna_se = mirna_se, class = "ER", pathway_num = 2)
names(table)
miRNA_path1 <- table[[1]]
Gene_path1 <- table[[2]]
```
The result would be a list, in matrix form, of related genes and miRNAs for each pathway. Note that because it would take a few minutes to run GSEA_ana, here we use the pre-calculated data to show the output.


5. GSEA_res

After doing GSEA analysis, we have selected miRNA and gene expression data for each enriched pathway. In the second step, the generated data would be put into “GSEA_res”.

```R
result <- GSEA_res(table, pheno.data = pheno.data, class = "ER", DE_method = "limma", cor_cut = -0.3)
names(result)
result_path1 <- result[[1]]
```

This function calculates the P-value, fold-change, and correlation for each miRNA-gene pair and shows these values to the user.


