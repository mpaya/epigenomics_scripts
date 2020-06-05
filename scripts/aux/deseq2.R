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

##################
## collect data ##
##################
fileData = read_tsv( "anno.dat", col_types = cols())

name_col1 <- names(fileData)[1]
n_order <- c(name_col1, "filenames", "sample")

sampleTable <- list.files(path = "./counts/", pattern = "*counts.txt") %>% 
    str_match("(.*).counts.txt") %>% as.data.frame %>% 
    dplyr::rename("filenames" = V1, !!name_col1 := V2) %>%
    left_join(fileData, by = setNames(name_col1, name_col1)) %>%
    select(c(n_order, setdiff(names(.), n_order))) 

dds <- DESeqDataSetFromHTSeqCount(sampleTable,
                                  directory = "./counts/",
                                  design = ~batch + sample )

cat("Total features:\t",nrow(dds),"\n", file="deseq.log", append=FALSE)
dds <- dds[ rowSums(counts(dds)) > 1, ]
cat("Features w/ counts:\t",nrow(dds),"\n\n", file="deseq.log", append=TRUE)

###########
## plots ##
###########
# Transform data with rlog
rld <- rlog(dds, blind = FALSE)

save(rld, dds, file = "rld.Rdata")
# load("rld.Rdata")

# plot with vst transformation
vsd <- vst(dds)
#pca_data <- plotPCA(vsd, intgroup = c("sample", "batch"), returnData=TRUE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
pca_data <- plotPCA(vsd, intgroup = c("sample", "batch"), returnData=TRUE)
p <- ggplot(pca_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = batch, shape = sample), size = 3) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1")
ggsave(p, file =  "pca_vsdBatch.tiff", width = 7, height = 7, compression="lzw")

# sample-sample distance matrix
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         filename = "sample_heatmap.tiff",
         width = 7, height = 7)

# PCA plot representation
pca_data <- plotPCA(rld, intgroup = c("sample", "batch"), returnData=TRUE)
p <- ggplot(pca_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = batch, shape = sample), size = 3) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1")
ggsave(p, file = "pca_rld.tiff", width = 7, height = 7)

#####################
## DESeq2 function ##
#####################
# dds <- DESeq(dds)
dds <- DESeq(dds, test = "LRT", reduced = ~batch)

write.table(counts( dds, normalized=TRUE ), file="deseq_normCnt.txt", sep = "\t", col.names=NA, quote = FALSE)
save(dds, file = "deseqAnalysis.RData")
# load(file = "deseqAnalysis.RData")

## collect results
ann <- read_tsv("genome_annot.txt")
sink("deseq.log", append=TRUE, split=TRUE) 
for (res_name in resultsNames(dds) %>% str_subset("sample")){
   s <- res_name %>% str_match("sample_(.*)_vs.*") %>% .[2]
   res <- results(dds, alpha=0.05 , name = res_name)
   cat(res_name,"\n")
   summary(res)
   res %>% as_tibble() %>% 
       add_column(genes = rownames(res), .before = 1) %>% 
       left_join(ann, by = c("genes" = "feature")) %>%
       write_tsv(paste0("deseqres_",s,"_annot.tsv"), na = "")
   cat(sapply(rownames(subset(res, padj < 0.05 & log2FoldChange<0)), toString), 
     file = paste0("DEtags_",s,"_down.txt", sep="\n"))
   cat(sapply(rownames(subset(res, padj < 0.05 & log2FoldChange>0)), toString), 
     file = paste0("DEtags_",s,"_up.txt", sep="\n"))
}
sink()

###################
### scatterplot ###
###################
## http://genoweb.toulouse.inra.fr/~pmartin/pgpmartin/2018/11/14/nicer-scatterplot-in-gggally/

if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse", repos='http://cran.us.r-project.org')
if (!requireNamespace("GGally", quietly = TRUE))
    install.packages("GGally", repos='http://cran.us.r-project.org')
if (!requireNamespace("fields", quietly = TRUE))
    install.packages("fields", repos='http://cran.us.r-project.org')

library(GGally)
library(tidyverse)

GGscatterPlot <- function(data, mapping, ..., 
                        method = "spearman") {

#Get correlation coefficient
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)

    cor <- cor(x, y, method = method)
#Assemble data frame
    df <- data.frame(x = x, y = y)
# PCA
    nonNull <- x!=0 & y!=0
    dfpc <- prcomp(~x+y, df[nonNull,])
    df$cols <- predict(dfpc, df)[,1]
# Define the direction of color range based on PC1 orientation:
    dfsum <- x+y
    colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                               dfsum[which.min(df$cols)],
                           1,
                           -1)
#Get 2D density for alpha
    dens2D <- MASS::kde2d(df$x, df$y)
    df$density <- fields::interp.surface(dens2D , 
                                         df[,c("x", "y")])

if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
}
#Prepare plot
    pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
                ggplot2::geom_point(shape=16, show.legend = FALSE) +
                ggplot2::scale_color_viridis_c(direction = colDirection) +
#                scale_color_gradient(low = "#0091ff", high = "#f0650e") +
                ggplot2::scale_alpha(range = c(.05, .6)) +
                ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
                ggplot2::geom_label(
                        data = data.frame(
                                        xlabel = min(x, na.rm = TRUE),
                                        ylabel = max(y, na.rm = TRUE),
                                        lab = round(cor, digits = 3)),
                        mapping = ggplot2::aes(x = xlabel, 
                                               y = ylabel, 
                                               label = lab),
                        hjust = 0, vjust = 1,
                        size = 3, fontface = "bold",
                        inherit.aes = FALSE # do not inherit anything from the ...
                        ) +
                theme_minimal()

return(pp)
}

df <- counts( dds, normalized=TRUE ) %>% as.data.frame

p <- GGally::ggpairs(log2(df+1),
                1:ncol(df),
                lower = list(continuous = GGscatterPlot),
                upper = list(continuous = wrap("cor", method= "spearman")))
ggsave(p, file="corr_rnaseq.tiff", width = ncol(df) + 1, height = ncol(df) + 1, compression="lzw")