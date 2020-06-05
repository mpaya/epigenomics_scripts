#!/usr/bin/Rscript

## install packages if needed and load
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
if (!requireNamespace("GenomicFeatures", quietly = TRUE)){
    BiocManager::install(c("GenomicRanges","GenomicFeatures"))}
if (!requireNamespace("ChIPpeakAnno", quietly = TRUE)){
    BiocManager::install("ChIPpeakAnno")}
if (!requireNamespace("tidyverse", quietly = TRUE)){
    install.packages("tidyverse")}
if (!requireNamespace("VennDiagram", quietly = TRUE)){
    install.packages("VennDiagram")}

suppressMessages(library(GenomicFeatures)); 
suppressMessages(library(ChIPpeakAnno)); 
suppressMessages(library(tidyverse));
suppressMessages(library(VennDiagram));

options(bitmapType='cairo')

##################
### annotation ###
##################

annotate_peaks <- function(bed_path, annoData){
    options(warn=-1)
    
    gr <- toGRanges(bed_path, format="broadPeak")

    peaks.anno <- annotatePeakInBatch(gr, AnnotationData=annoData, 
                                  output="overlapping")
    peak_table <- as.data.frame(peaks.anno)
    
    options(warn=0)
    return(peak_table)
}


path = getwd()

pattern <- "epic2_res_(.*)_(\\d+).bed"
res_tbl <- list.files(path, pattern = "*\\.bed")  %>% 
    str_match(pattern) %>% as_tibble() %>% na.omit %>%
    rename("beds" = 1, "sample" = 2, "gap" = 3)
## attempt to get total peaks from raw res
pattern2 <- "epic2_res_(.*)_(\\d+).txt"
res_tbl <- list.files(path, pattern = "*\\.txt")  %>% 
    str_match(pattern2) %>% as_tibble() %>% na.omit %>%
    rename("rawres" = 1, "sample" = 2, "gap" = 3) %>%
    full_join(res_tbl, by = c("sample","gap"))


gff <- list.files(path, pattern = "\\.gff") 
txdb <- makeTxDbFromGFF(file=gff, format="gff", dataSource="Araport11", organism="Arabidopsis")
annoData <- toGRanges(txdb, feature="gene")


gene_lists <- list()
df <- data.frame(cond = character(1), gap = character(1), n_islands = numeric(1), 
	             filt_pk = numeric(1), mean_len = numeric(1), 
	             median = numeric(1), min_len = numeric(1), max_len = numeric(1), 
				 n_genes = character(1), score = numeric(1), stringsAsFactors=FALSE)
largest_peaks <- data.frame()

for (idx in seq(nrow(res_tbl))) {
	bed_file <- res_tbl$beds[idx]
	raw_peaks <- res_tbl$rawres[idx] %>% read_tsv(col_types = cols()) %>% nrow
	samp <- paste(res_tbl$sample[idx], res_tbl$gap[idx], sep='-')
	annot_res <- annotate_peaks(bed_file, annoData)
	
	gene_lists[[samp]] <- pull(annot_res, feature) %>% unique %>% na.omit
	num_genes <- gene_lists[[samp]] %>% length
	
	raw_epic <- annot_res %>% select(c(1:10)) %>% unique
	
	df <- add_row(df, cond = res_tbl$sample[idx], gap = res_tbl$gap[idx], 
		          n_islands = raw_peaks, filt_pk = nrow(raw_epic), 
				  mean_len = mean(raw_epic$width), median = median(raw_epic$width), 
				  min_len = min(raw_epic$width), max_len = max(raw_epic$width),
				  n_genes = num_genes, score = sum(raw_epic$score))
	
	largest_peaks <- raw_epic %>% arrange(desc(width)) %>% 
			add_column(file = bed_file, .before = 1) %>% 
			head(10) %>% bind_rows(largest_peaks, .)
}
df <- df[-1,]

print("results collected")

# Plot aggregate score as a function of the gap size

p <- ggplot() + 
   geom_text(data=df, aes(x = mean_len, y = score, color = cond, label = gap), size=5) + 
   scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
   labs(x= "Mean length", y = "Aggregated score")

ggsave(p, file = "epic_gapSize.tiff", width = 5, height = 4, compression="lzw")

print("aggregate score plotted")

# Add coverage res
cov <- list.files('.', pattern = "*_coverage.tsv")
cov <- read_tsv(cov, col_names = FALSE, col_types = cols(), comment = "#") %>% 
		select(c(1,cond = 2, gap = 3,6)) %>% spread(1,4) %>% 
		mutate(gap = as.character(gap))
df <- left_join(df, cov, by = c("cond", "gap"))

## Jaccard of results
jac_res <- list.files('.', pattern = "*_jaccard.tsv")
jac_res <- read_tsv(jac_res, col_types = cols())
names(jac_res) <- gsub("[-_]",".", names(jac_res))
p <- ggplot(data=jac_res, aes_string(x = names(jac_res)[1], y = names(jac_res)[2])) + 
    geom_tile(aes(fill = jaccard)) + 
    geom_text(aes(label = round(jaccard, 3)))
ggsave(p, file = "epic_gapJaccard.tiff", width=6, height=4, compression="lzw")

print("jaccard done")

## Venn diagrams
s_list <- df %>% pull(cond) %>% unique

for (i in 1:length(s_list)) {
	venn.plot <- VennDiagram::venn.diagram(gene_lists[paste(
            s_list[i], unique(df$gap)[1:5], sep = "-")], 
    cex = 1, cat.cex = 1.5, margin = .1, filename = paste0("venn_",s_list[i],".tiff"))
}
print("venn done")

# Save results
df <- sapply(df, as.character)
write.table(df, "epic_gapRes.tsv", quote=FALSE, 
            row.names=FALSE, sep="\t", na = "")
write_tsv(largest_peaks, "largest_peaks.txt")

print("process finished")
