


library(tidyverse)


sample_names <- read.table(file.path("D:/2_POST_DOC/HKO_SEQ_DATA/HKO13045-98258"
                                     , "sample_names.txt"), header = T)

sample_names %>% head


dir <- "D:/2_POST_DOC/HKO_SEQ_DATA/HKO13045-98258/Kallisto"

files <- file.path(dir, sample_names$samples, "abundance.tsv")

names(files) <- c("EMT6", "EMT6-DRG", "4T1", "4T1-DRG")

files %>% print

### I will prepare a transcript to gene data frame (tx2gene)

library(ensembldb)

mouse_ens_v108 <- 
  makeTxDbFromGFF("D:/2_POST_DOC/HKO_SEQ_DATA/Mus_musculus.GRCm39.108.gtf",
                  format="gtf")

### I will look at all the keytypes available:

keytypes(mouse_ens_v108)

### I will view the some keytypes of this database object, by keys():
 
keys(mouse_ens_v108, keytype = "GENEID") %>% head()

keys(mouse_ens_v108, keytype = "GENEID") %>% length()

### It contains the ensembl ids for 56980 genes.

keys(mouse_ens_v108, keytype = "TXNAME") %>% head()

keys(mouse_ens_v108, keytype = "TXNAME") %>% length()

### It contains 149482 ensemble transcript ids.

### I will retrieve ensmebl gene id to transcript id conversion data:

df_mouse_ens_v108 <- select(mouse_ens_v108, keys = keys(mouse_ens_v108, 
                                                    keytype = "GENEID"), 
                                                    keytype = "GENEID", 
                                                    columns = "TXNAME")


df_mouse_ens_v108 %>% head

df_mouse_ens_v108 %>% class()

tx2gene_mouse_ens_v108 <- df_mouse_ens_v108[, 2:1]  # tx ID, then gene ID

tx2gene_mouse_ens_v108 %>% head

### Kallisto .tsv data includes the additional number at the end fo ensemble
### transcript id that displays the transcript version, separated by ".".
### We ignored this by setting the "ignoreTxVersion" as "TRUE".

txi.kallisto.tsv <- tximport::tximport(files, type = "kallisto", 
                             tx2gene = tx2gene_mouse_ens_v108,
                             ignoreTxVersion = TRUE)

### Exploratory Analysis:

### First, we prepare a sample table to supply to deseq2 function:

sample_table <- data.frame(sample_names=colnames(txi.kallisto.tsv$abundance),
                condition=factor(x = c("control", "wDRG", "control","wDRG")),
                           cell_line=factor(x=c("EMT6","EMT6","4T1","4T1")))

### Our conditions are categorical variables so they should be factors.

sample_table %>% print

library(DESeq2)

dds_cancer_DRG <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                        sample_table, ~ cell_line + condition)

### In the multifactor design,  As condition is the variable of interest, 
### we put it at the end of the formula. We can access to the design formula
### and change it by design function:

design(dds_cancer_DRG) %>% print()

dds_cancer_DRG %>% nrow()

### Filtering the rows less than 10 in sum:

keep <- rowSums(counts(dds_cancer_DRG)) >= 10

dds_cancer_DRG <- dds_cancer_DRG[keep, ]

nrow(dds_cancer_DRG)

### Regularized log2 transformation:

rld_b_f <- rlog(dds_cancer_DRG, blind=FALSE)

plotPCA(rld_b_f, intgroup=c("condition"))

plotPCA(rld_b_f, intgroup=c("sample_names"))

plotPCA(rld_b_f, intgroup=c("cell_line"))

### Sample distances. Default function is euclidean.

sampleDists <- dist(t(assay(rld_b_f)))

library("pheatmap")
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(as.matrix(sampleDists), 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
                                   col = colors)

### one can make a boxplot of the Cook’s distances to see if one sample is 
### consistently higher than others (here this is not the case):

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_cancer_DRG)[["cooks"]]), range=0, las=2)

# Plotting the dispersion estimates is a useful diagnostic. The dispersion 
# plot below is typical, with the final estimates shrunk from the gene-wise 
# estimates towards the fitted estimates. Some gene-wise estimates are flagged 
# as outliers and not shrunk towards the fitted value, (this outlier detection 
# is described in the manual page for estimateDispersionsMAP). The amount of 
# shrinkage can be more or less than seen here, depending on the sample size, 
# the number of coefficients, the row mean and the variability of the gene-wise 
# estimates.

plotDispEsts(dds_cancer_DRG)

### Differential Expression Analysis:

dds_cancer_DRG <- DESeq(dds_cancer_DRG)

res <- 
  results(dds_cancer_DRG, 
          contrast = c("condition", "wDRG", "control"), alpha = 0.05)



### how many significant genes we have:

sum(res$padj < 0.05, na.rm=TRUE)


### MA plot:

plotMA(res, ylim=c(-8,8))

### add names to the all results table:

library(org.Mm.eg.db)

res$symbols <- mapIds(org.Mm.eg.db, keys = row.names(res), 
                          keytype = "ENSEMBL", column = "SYMBOL",
                          multiVals = "first")

res$description <- mapIds(org.Mm.eg.db, keys = row.names(res), 
                              keytype = "ENSEMBL", column = "GENENAME",
                              multiVals = "first")

### order by padj:

res <- res[order(res$padj, decreasing = FALSE),]

### significant genes:

res_sig <- subset(x = res, res$padj < 0.05)

### filling  two rows' gene names and description:

res_sig$symbols["ENSMUSG00000083365"] <- "Gm12895-201"
res_sig$symbols["ENSMUSG00000083558"] <- "Gm12896-201"

res_sig$description["ENSMUSG00000083365"]  <- "Processed pseudogene"
res_sig$description["ENSMUSG00000083558"]  <- "Processed pseudogene"

### info about the properties:

mcols(res, use.names=TRUE)

### The first column, baseMean, is a just the average of the normalized count 
### values, dividing by size factors, taken over all samples in the 
### DESeqDataSet.

### heatmap:

mat <- assay(rld_b_f)[ row.names(res_sig), ]
df <- as.data.frame(colData(rld_b_f)[,c("condition","cell_line")])

### I will not use factor (1,2) in the colData but the cell line names
### in the heatmap.


### prepare a dataframe of gene names to use its column as row labels:

annot_row_df <- res_sig$symbols[row.names(mat)] %>% as.data.frame()

pheatmap::pheatmap(mat, scale = "row", annotation_col=df, 
                   labels_row = annot_row_df[,1])


### Alternative heatmap function:

# gplots::heatmap.2(mat, 
#           scale = "row", trace="none", 
#           col = gplots::colorpanel(64,"red", "black", "green"), 
#           labRow = annot_row_df[,1] ,
#           margins = c(16, 21) ,dendrogram = "none",  
#           cexCol = 0.5, cexRow = 0.5)

### Prepare the FPKM values:

# The following function returns fragment counts normalized per kilobase of 
# feature length per million mapped fragments (by default using a robust 
# estimate of the library size, as in estimateSizeFactors).


fpkm_dds <- fpkm(object = dds_cancer_DRG, robust = TRUE)

### default:  robust = TRUE 

fpkm_dds <- as.data.frame(fpkm_dds)

fpkm_dds %>% dim()

res %>% dim()

### I will add the gene name, padj and description columns by row.names:

fpkm_w_padj <- merge(x = as.data.frame(res), y = fpkm_dds, by="row.names")

fpkm_w_padj %>% head

fpkm_w_padj[,c(2,4:6)] <- NULL

### change the order of the columns:

fpkm_w_padj <- fpkm_w_padj[,c(1,4:9,3,2)]

### order by adjusted p value:

fpkm_w_padj <- fpkm_w_padj[order(fpkm_w_padj$padj, decreasing = FALSE),]

### export:

write.csv(x = res_sig, file = "significant_genes.csv")
write.csv(x = res, file = "results.csv")

write.csv(x = fpkm_w_padj, file = "fpkm_padj.csv")

### 

### Normalized counts are not log2 transformed.

normalised_counts <- counts(dds_cancer_DRG, normalized = TRUE)

normalised_counts$symbols <-  mapIds(org.Mm.eg.db, 
                                     keys = row.names(normalised_counts), 
                                     keytype = "ENSEMBL", column = "SYMBOL",
                                     multiVals = "first") 

normalised_counts$description <-  mapIds(org.Mm.eg.db, 
                                         keys = row.names(normalised_counts), 
                                         keytype = "ENSEMBL", column = "GENENAME",
                                         multiVals = "first") 

write.csv(x = normalised_counts, file = "normalised_counts.txt")

### After retrieval of the normalized counts from the deseq2 object,
### I calculated log2fold change values manually for comparisons
### of samples, each having one replicate. I did this manually because
### new version of the deseq2 package does not allow the comparison with
### only single sample in each group.

### Log2(A/B) <- Log2(A) - Log2(B). So I used this equation for the
### following calculations.

normalised_counts$'4T1_vs_EMT6' <- 
  log2(1+(normalised_counts$`4T1`)) - log2(1+(normalised_counts$EMT6))

normalised_counts$'4T1_DRG_vs_4T1' <- 
  log2(1+(normalised_counts$`4T1-DRG`)) - log2(1+(normalised_counts$`4T1`))

normalised_counts$'EMT6_DRG_vs_EMT6' <- 
  log2(1+(normalised_counts$`EMT6-DRG`)) - log2(1+(normalised_counts$EMT6))

### Form fpkm data:

fpkm_dds <- fpkm(object = dds_cancer_DRG, robust = TRUE)
fpkm_dds <- as.data.frame(fpkm_dds)
fpkm_dds %>% dim()
 
fpkm_w_LogFC_no_replicate <- merge(x = fpkm_dds, 
                                   y = normalised_counts[,-c(1:4)], 
                                   by="row.names")

### I exclude the first four column, since they have the same column names
### with the fpkm df column names, and also we do not need the counts data.

fpkm_w_LogFC_no_replicate %>% head()

###

### gene list for log2FC larger than 2, so UP in 4T1 vs EMT6:

fpkm_4T1_vs_EMT6 <- subset(fpkm_w_LogFC_no_replicate, 
                           `4T1_vs_EMT6` > 2 )

### I subset to those with fpkm > 0.5 in 4T1 sample:

fpkm_4T1_vs_EMT6 <- subset(fpkm_4T1_vs_EMT6, `4T1` > 0.5)

### gene list for log2FC smaller than -2, so UP in EMT6 vs 4T1:

fpkm_EMT6_vs_4T1 <- subset(fpkm_w_LogFC_no_replicate, 
                           `4T1_vs_EMT6` < -2 )

fpkm_EMT6_vs_4T1 <- subset(fpkm_EMT6_vs_4T1, EMT6 > 0.5)

### One should not be confused by my usage of < -2 . I manually calculated
### the log2 fold changes only for one direction, A/B but not B/A. Therefore,
### To get B/A, I use the minus sign.

###

fpkm_4T1_DRG_vs_4T1 <- subset(fpkm_w_LogFC_no_replicate, 
                      `4T1_DRG_vs_4T1` > 2 )

fpkm_4T1_DRG_vs_4T1 <- subset(fpkm_4T1_DRG_vs_4T1, `4T1-DRG` > 0.5)

fpkm_4T1_DRG_vs_4T1 %>% dim()

###

fpkm_EMT6_DRG_vs_EMT6 <- subset(fpkm_w_LogFC_no_replicate, 
                                EMT6_DRG_vs_EMT6 > 2 )

fpkm_EMT6_DRG_vs_EMT6 <- subset(fpkm_EMT6_DRG_vs_EMT6, 
                                `EMT6-DRG` > 0.5)

fpkm_EMT6_DRG_vs_EMT6 %>% dim()

###

fpkm_4T1_vs_4T1_DRG <- subset(fpkm_w_LogFC_no_replicate, 
                              `4T1_DRG_vs_4T1` < -2 )

fpkm_4T1_vs_4T1_DRG <- subset(fpkm_4T1_vs_4T1_DRG, `4T1` > 0.5)

fpkm_4T1_vs_4T1_DRG %>% dim()

###

### I also placed another threshold for the fpkm values otherwise 
### log2 fold changes might be misleading especially for the lowly expressed
### genes.

###

fpkm_EMT6_vs_EMT6_DRG <- subset(fpkm_w_LogFC_no_replicate, 
                                EMT6_DRG_vs_EMT6 < -2 )

fpkm_EMT6_vs_EMT6_DRG <- subset(fpkm_EMT6_vs_EMT6_DRG, 
                                EMT6 > 0.5)

fpkm_EMT6_vs_EMT6_DRG %>% dim()

### I will add the gene symbols and descriptions using R interface of gprofiler:

###The previous method that I used returned so many NAs.

### NOTE: Previously, I merged two data frames by row.names , which created a
### a first column with row.names. I did not change this first column named
### named as "row.names". So, I refer this column of the datbase 
### with "$row.names".

fpkm_4T1_DRG_vs_4T1$symbols <- 
  gprofiler2::gconvert(query = fpkm_4T1_DRG_vs_4T1$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

fpkm_4T1_DRG_vs_4T1$description <- 
  gprofiler2::gconvert(query = fpkm_4T1_DRG_vs_4T1$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$description

###

fpkm_4T1_vs_4T1_DRG$symbols <- 
  gprofiler2::gconvert(query = fpkm_4T1_vs_4T1_DRG$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

fpkm_4T1_vs_4T1_DRG$description <- 
  gprofiler2::gconvert(query = fpkm_4T1_vs_4T1_DRG$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$description

###

fpkm_4T1_vs_EMT6$symbols <- 
  gprofiler2::gconvert(query = fpkm_4T1_vs_EMT6$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

fpkm_4T1_vs_EMT6$description <- 
  gprofiler2::gconvert(query = fpkm_4T1_vs_EMT6$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$description




###

fpkm_EMT6_DRG_vs_EMT6$symbols <- 
  gprofiler2::gconvert(query = fpkm_EMT6_DRG_vs_EMT6$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

fpkm_EMT6_DRG_vs_EMT6$description <- 
  gprofiler2::gconvert(query = fpkm_EMT6_DRG_vs_EMT6$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$description

###

fpkm_EMT6_vs_4T1$symbols <- 
  gprofiler2::gconvert(query = fpkm_EMT6_vs_4T1$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

### It returned error:

### Error in `$<-.data.frame`(`*tmp*`, symbols, value = c("Axin2", "Gmpr",: 

### replacement has 820 rows, data has 821.

### To see what is missing in the larger set:

base::setdiff(fpkm_EMT6_vs_4T1$Row.names, my_df$input)

#### "ENSMUSG00000121513"

### From ensembl website, symbol: PTPRG. 
### Description:  protein tyrosine phosphatase, receptor type, G.

my_df <- 
  gprofiler2::gconvert(query = fpkm_EMT6_vs_4T1$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG")

fpkm_EMT6_vs_4T1$symbols <- c(my_df$name, "PTPRG") ## I also added the last one.

fpkm_EMT6_vs_4T1$symbols %>% tail()

fpkm_EMT6_vs_4T1$description <- 
  c(my_df$description,  "protein tyrosine phosphatase, receptor type, G") 
 
fpkm_EMT6_vs_4T1$description %>% tail()

fpkm_EMT6_vs_4T1$symbols %>% head
### I also added the last one.

fpkm_EMT6_vs_EMT6_DRG$symbols <- 
  gprofiler2::gconvert(query = fpkm_EMT6_vs_EMT6_DRG$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$name

fpkm_EMT6_vs_EMT6_DRG$description <- 
  gprofiler2::gconvert(query = fpkm_EMT6_vs_EMT6_DRG$Row.names,
                       organism = "mmusculus", 
                       target = "ENSG") %>% .$description

fpkm_EMT6_vs_EMT6_DRG$description %>% length()

fpkm_EMT6_vs_EMT6_DRG$description %>% print

###


###
write.csv(fpkm_4T1_DRG_vs_4T1, file = "fpkm_4T1_DRG_vs_4T1.csv")

write.csv(fpkm_4T1_vs_4T1_DRG, file = "fpkm_4T1_vs_4T1_DRG.csv")

write.csv(fpkm_4T1_vs_EMT6, file = "fpkm_4T1_vs_EMT6.csv")

write.csv(fpkm_EMT6_DRG_vs_EMT6, file = "fpkm_EMT6_DRG_vs_EMT6.csv")

write.csv(fpkm_EMT6_vs_4T1, file = "fpkm_EMT6_vs_4T1.csv")

write.csv(fpkm_EMT6_vs_EMT6_DRG, file = "fpkm_EMT6_vs_EMT6_DRG.csv")


###

mat_all <- assay(rld_b_f)

mat_all %>% colnames()
### "EMT6"     "EMT6-DRG" "4T1"      "4T1-DRG" 


### heatmaps:

labels_row <- data.frame(genes=fpkm_EMT6_vs_EMT6_DRG$symbols,
                        row.names = fpkm_EMT6_vs_EMT6_DRG$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_EMT6_vs_EMT6_DRG$Row.names,], 
          scale = "row", 
          annotation_col=df, 
          labels_row = labels_row)

### Heatmap with only EMT and EMT-DRG columns:


pheatmap::pheatmap(mat_all[fpkm_EMT6_vs_EMT6_DRG$Row.names,1:2],
                   scale = "row",
                   annotation_col=df,
                   labels_row = labels_row)

### It did not return a meaningful view. Unscaling rows can also be tried:

# pheatmap::pheatmap(mat_all[fpkm_EMT6_vs_EMT6_DRG$Row.names,1:2], 
#                    scale = "none", 
#                    annotation_col=df, 
#                    labels_row = labels_row)

### Again, not interpretable.

### Heatmaps for all other pair wise comparisons:

labels_row <- data.frame(genes=fpkm_4T1_DRG_vs_4T1$symbols,
                         row.names = fpkm_4T1_DRG_vs_4T1$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_4T1_DRG_vs_4T1$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 3)


labels_row <- data.frame(genes=fpkm_4T1_vs_4T1_DRG$symbols,
                         row.names = fpkm_4T1_vs_4T1_DRG$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_4T1_vs_4T1_DRG$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 3)

labels_row <- data.frame(genes=fpkm_4T1_vs_EMT6$symbols,
                         row.names = fpkm_4T1_vs_EMT6$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_4T1_vs_EMT6$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 3)

labels_row <- data.frame(genes=fpkm_EMT6_vs_4T1$symbols,
                         row.names = fpkm_EMT6_vs_4T1$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_EMT6_vs_4T1$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 3)

labels_row <- data.frame(genes=fpkm_EMT6_vs_EMT6_DRG$symbols,
                         row.names = fpkm_EMT6_vs_EMT6_DRG$Row.names) %>% .[,1]

pheatmap::pheatmap(mat_all[fpkm_EMT6_vs_EMT6_DRG$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 8)

labels_row <- data.frame(genes=fpkm_EMT6_DRG_vs_EMT6$symbols,
                         row.names = fpkm_EMT6_DRG_vs_EMT6$Row.names) %>% .[,1]


pheatmap::pheatmap(mat_all[fpkm_EMT6_DRG_vs_EMT6$Row.names,], 
                   scale = "row", 
                   annotation_col=df, 
                   labels_row = labels_row, fontsize_row = 5)

sessionInfo()

# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils   datasets  methods   base     
# 
# other attached packages:
#   [1] org.Mm.eg.db_3.15.0         DESeq2_1.36.0              
# [3] SummarizedExperiment_1.26.1 MatrixGenerics_1.8.1       
# [5] matrixStats_0.62.0          ensembldb_2.20.2           
# [7] AnnotationFilter_1.20.0     GenomicFeatures_1.48.4     
# [9] AnnotationDbi_1.58.0        Biobase_2.56.0             
# [11] GenomicRanges_1.48.0        GenomeInfoDb_1.32.4        
# [13] IRanges_2.30.1              S4Vectors_0.34.0           
# [15] BiocGenerics_0.42.0         forcats_0.5.2              
# [17] stringr_1.4.1               dplyr_1.0.10               
# [19] purrr_0.3.5                 readr_2.1.3                
# [21] tidyr_1.2.1                 tibble_3.1.8               
# [23] ggplot2_3.3.6               tidyverse_1.3.2            
# 
# loaded via a namespace (and not attached):
#   [1] googledrive_2.0.0        colorspace_2.0-3         rjson_0.2.21            
# [4] ellipsis_0.3.2           XVector_0.36.0           fs_1.5.2                
# [7] rstudioapi_0.14          farver_2.1.1             bit64_4.0.5             
# [10] fansi_1.0.3              lubridate_1.8.0          xml2_1.3.3              
# [13] splines_4.2.1            codetools_0.2-18         tximport_1.24.0         
# [16] cachem_1.0.6             geneplotter_1.74.0       jsonlite_1.8.2          
# [19] Rsamtools_2.12.0         broom_1.0.1              annotate_1.74.0         
# [22] dbplyr_2.2.1             png_0.1-7                pheatmap_1.0.12         
# [25] compiler_4.2.1           httr_1.4.4               backports_1.4.1         
# [28] assertthat_0.2.1         Matrix_1.5-1             fastmap_1.1.0           
# [31] lazyeval_0.2.2           gargle_1.2.1             cli_3.4.1               
# [34] prettyunits_1.1.1        tools_4.2.1              gtable_0.3.1            
# [37] glue_1.6.2               GenomeInfoDbData_1.2.8   rappdirs_0.3.3          
# [40] Rcpp_1.0.9               cellranger_1.1.0         vctrs_0.4.2             
# [43] Biostrings_2.64.1        rhdf5filters_1.8.0       rtracklayer_1.56.1      
# [46] rvest_1.0.3              lifecycle_1.0.3          restfulr_0.0.15         
# [49] XML_3.99-0.11            googlesheets4_1.0.1      zlibbioc_1.42.0         
# [52] scales_1.2.1             vroom_1.6.0              hms_1.1.2               
# [55] ProtGenerics_1.28.0      parallel_4.2.1           rhdf5_2.40.0            
# [58] RColorBrewer_1.1-3       yaml_2.3.5               curl_4.3.3              
# [61] memoise_2.0.1            biomaRt_2.52.0           stringi_1.7.8           
# [64] RSQLite_2.2.18           genefilter_1.78.0        BiocIO_1.6.0            
# [67] filelock_1.0.2           BiocParallel_1.30.4      rlang_1.0.6             
# [70] pkgconfig_2.0.3          bitops_1.0-7             lattice_0.20-45         
# [73] Rhdf5lib_1.18.2          labeling_0.4.2          GenomicAlignments_1.32.1
# [76] bit_4.0.4                tidyselect_1.2.0         magrittr_2.0.3          
# [79] R6_2.5.1                 generics_0.1.3           DelayedArray_0.22.0     
# [82] DBI_1.1.3                pillar_1.8.1             haven_2.5.1             
# [85] withr_2.5.0              survival_3.4-0           KEGGREST_1.36.3         
# [88] RCurl_1.98-1.9           modelr_0.1.9             crayon_1.5.2            
# [91] utf8_1.2.2               BiocFileCache_2.4.0      tzdb_0.3.0              
# [94] progress_1.2.2           locfit_1.5-9.6           grid_4.2.1              
# [97] readxl_1.4.1             blob_1.2.3               reprex_2.0.2            
# [100] digest_0.6.29            xtable_1.8-4             munsell_0.5.0  