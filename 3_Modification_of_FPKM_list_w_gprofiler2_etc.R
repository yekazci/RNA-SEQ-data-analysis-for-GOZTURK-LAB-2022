


setwd("D:/2_POST_DOC/HKO_SEQ_DATA/HKO13045-98258/Kallisto/s1_1_2_IDT8_UDI_179")
load("D:/2_POST_DOC/HKO_SEQ_DATA/HKO13045-98258/Kallisto/Kallisto_tximport_logFC_w_no_replicates.RData")
library(tidyverse)
library(DESeq2)
library(gprofiler2)

fpkm_dds %>% dim()
### [1] 18123     4

fpkm_w_symbols_etc <- fpkm_dds

my_df$input %>% length()

### One ensembl id returned no gene symbol. gProfiler2 deletes the input id
### if there is not any symbol.

### I will use setdiff from the R base and identify the ensembl id in the 
### gprofiler result df that differs from the input data:
### 

base::setdiff(row.names(fpkm_w_symbols_etc), my_df$input)

#### "ENSMUSG00000121513"

### From ensembl website, symbol: PTPRG. 
### Description:  protein tyrosine phosphatase, receptor type, G.


base::setdiff(my_df$input, row.names(fpkm_w_symbols_etc))

### character(0)

### select the columns of interest:

my_df <- my_df[,c("input", "name", "description")]

### adding the missing ensembl id to the end row:

my_df[(nrow(my_df)+1), 1:3] <- c("ENSMUSG00000121513", "PTPRG", 
                              "protein tyrosine phosphatase, receptor type, G")

my_df %>% tail()

fpkm_w_symbols_etc %>% head()

### Before merging, I will set the row names as a new column.

fpkm_w_symbols_etc$input <- row.names(fpkm_w_symbols_etc)

### Now, both data frames have ensembl ids in the column named "input"

### I will merge them by the "input".

fpkm_w_symbols_etc <- merge(fpkm_w_symbols_etc, my_df, by="input")

fpkm_w_symbols_etc %>% head()

fpkm_w_symbols_etc$name %>% anyNA()

### FALSE

res_temp <- res

res_temp$input <- row.names(res_temp)

res_temp <- as.data.frame(res_temp[,c("input", "log2FoldChange")])

### merge by :

fpkm_w_symbols_etc <- merge(fpkm_w_symbols_etc, 
                            res_temp[,c("input", "log2FoldChange")], by="input")

fpkm_w_symbols_etc %>% head()

### change the first column name into a more sensible name:

colnames(fpkm_w_symbols_etc)[1] <- "ensembl_ids"

fpkm_w_symbols_etc %>% head()

remove(res_temp)

library("SBGNview")

xx <- data("cancer.ds")

remove(xx)

fpkm_w_symbols_etc %>% head()

my_data <- fpkm_w_symbols_etc[, "ensembl_ids", "log2FoldChange"]

### pathview uses either a named vector or data.frame/matrix with row.names.

### names or row.names must be gene ids and column or vector elements must be
### log 2 fold changes.

pth_obj <- pathview(gene.data = my_data, pathway.id = "04657",species = "mmu", 
         out.suffix = "cancer_DRG",gene.idtype = "ensembl",
         low = list(gene = "red", cpd = "blue"), 
         mid = list(gene = "gray", cpd = "gray"), 
         high = list(gene = "green", cpd ="yellow"))

pth_obj$plot.data.gene %>% head() ### we can access the pathway mapped genes.

write.csv(pth_obj$plot.data.gene, file = "mmu04657_pathview_results.csv")

###

write.csv(fpkm_w_symbols_etc, file = "fpkm_w_symbols_etc.csv")
