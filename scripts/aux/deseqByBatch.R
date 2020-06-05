#!/usr/bin/Rscript

if (!requireNamespace("tidyverse", quietly = TRUE)){
    install.packages("tidyverse")}
if (!requireNamespace("pheatmap", quietly = TRUE)){
    install.packages("pheatmap")}
if (!requireNamespace("RColorBrewer", quietly = TRUE)){
    install.packages("RColorBrewer")}
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
if (!requireNamespace("tidyverse", quietly = TRUE)){
    BiocManager::install("DESeq2")}

suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))

options(bitmapType='cairo')

## collect data
fileData = read_tsv( "anno.dat", col_types = "ffccccf")
ann <- read_tsv("genome_annot.txt")

name_col1 <- names(fileData)[1]
n_order <- c(name_col1, "filenames", "sample")
sampleTable <- list.files(path = "./counts/", pattern = "*counts.txt") %>% 
    str_match("(.*).counts.txt") %>% as.data.frame %>% 
    dplyr::rename("filenames" = V1, !!name_col1 := V2) %>%
    left_join(fileData, by = setNames(name_col1, name_col1)) %>%
    select(c(n_order, setdiff(names(.), n_order))) 

cat("DESeq2 by batch. \n", file="batch_deseq.log", append=FALSE)

for (b in unique(sampleTable$batch)){
  sT <- filter(sampleTable, batch == b)
  sT$sample <- factor(sT$sample, levels = unique(sT$sample)) %>% relevel(ref = "WT")

  dds <- DESeqDataSetFromHTSeqCount(sT,
                                  directory = "./counts/",
                                  design = ~ sample )

  cat(b," batch. Total features:\t",nrow(dds),"\n", file="batch_deseq.log", append=T)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  cat(b," batch. Feat w/ counts:\t",nrow(dds),"\n\n", file="batch_deseq.log", append=T)


  ## DESeq2 function
  dds <- DESeq(dds, test = "LRT", reduced = ~1)

  ## collect results by sample
  for (res_name in resultsNames(dds)[2:length(resultsNames(dds))]) {
    s <- str_split(res_name, "_")[[1]][2]
    res <- results(dds, alpha=0.05 , name = res_name)
    sink("batch_deseq.log", append=TRUE, split=TRUE)
    cat(res_name); summary(res)
    sink()
    res %>% as_tibble() %>% 
        add_column(genes = rownames(res), .before = 1) %>% 
        left_join(ann, by = c("genes" = "feature")) %>%
        write_tsv(paste0("bybatchRes_",s,"_annot.tsv"), na = "")
  }
}

