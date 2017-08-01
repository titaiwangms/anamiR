# anamiR-2.0
An integrated analysis R package of miRNA and mRNA profiiling
## anamiR
### General Workflow
1. Import data
```R
library(anamiR)
data(mrna)
data(mirna)
data(pheno.mirna)
data(pheno.mrna)
```
2. SummarizedExperimaent object
```R
mirna_se <- SummarizedExperiment(assays = SimpleList(counts=mirna), colData = pheno.mirna)
mrna_se <- SummarizedExperiment(assays = SimpleList(counts=mrna), colData = pheno.mrna)
```
3. Differential Expression analysis
```R
mirna_d <- differExp_discrete(se = mirna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
mrna_d <- differExp_discrete(se = mrna_se, class = "ER", method = "limma", log2 = FALSE, p_value.cutoff = 0.05, logratio = 0.5, p_adjust.method = "BH")
```
4. miRNA ID conversion
```R
mirna_21 <- miR_converter(data = mirna_d, remove_old = TRUE, original_version = 17, latest_version = 21)
```
5. Correlation abalysis
```R
cor <- negative_cor(mrna_data = mrna_d, mirna_data = mirna_21, cut.off = -0.5, method = "pearson")
```
6. Database Intersection
```R
sup <- database_support(cor_data = cor, Sum.cutoff = 2, org = "hsa")
```
7. Heatmap visualization
```R
heat_vis(cor, mrna_d, mirna_21)
```
8. Functional analysis
```R
enr <- enrichment(data_support = sup, per_time = 5000)
```
### GSEA Workflow
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
```R
table <- GSEA_ana(mrna_se = mrna_se, mirna_se = mirna_se, class = "ER", pathway_num = 2)
names(table)
miRNA_path1 <- table[[1]]
Gene_path1 <- table[[2]]
```
5. GSEA_res
```R
result <- GSEA_res(table, pheno.data = pheno.data, class = "ER", DE_method = "limma", cor_cut = -0.3)
names(result)
result_path1 <- result[[1]]
```
