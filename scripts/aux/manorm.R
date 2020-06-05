#!/usr/bin/Rscript
args <- commandArgs(trailingOnly=TRUE)

## install packages if needed and load
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
if (!requireNamespace("GenomicFeatures", quietly = TRUE)){
    BiocManager::install("GenomicFeatures")}
if (!requireNamespace("ChIPpeakAnno", quietly = TRUE)){
    BiocManager::install("ChIPpeakAnno")}
if (!requireNamespace("tidyverse", quietly = TRUE)){
    install.packages("tidyverse")}

suppressMessages(library(GenomicFeatures)); 
suppressMessages(library(ChIPpeakAnno)); 
suppressMessages(library(tidyverse));

## functions
get_metadat <- function(bed, gr1, name = "s"){
    gr2 <- toGRanges(bed, format="broadPeak")
    ol <- findOverlaps(gr1, gr2)

    peaks <- cbind(as_tibble(ol), peaks = names(gr2)[as.matrix(ol)[,2]])
    peaks <- aggregate(peaks ~ queryHits, data = peaks, paste, collapse = ",")
    peaks <- left_join(tibble(queryHits = seq(nrow(values(gr1)))), peaks, by = "queryHits")
    
    FC <- cbind(as_tibble(ol), fc = gr2$signalValue[as.matrix(ol)[,2]])
    FC <- aggregate(fc ~ queryHits, data = FC, paste, collapse = ",")
    FC <- left_join(tibble(queryHits = seq(nrow(values(gr1)))), FC, by = "queryHits")
    
    meta_dat <- tibble(!!paste0("peaks_",name) := peaks$peaks, 
                       !!paste0("FC_",name) := FC$fc) 
    return(meta_dat)
}

annotate_peaks <- function(gr, annoData){
    options(warn=-1)
    
    if (is.character(gr)){
        gr <- toGRanges(gr, format="broadPeak")
    }

    peaks.anno <- annotatePeakInBatch(gr, AnnotationData=annoData, 
                                  output="overlapping", maxgap=500L)
    peak_table <- as.data.frame(peaks.anno)
    
    options(warn=0)
    return(peak_table)
}

write_annot <- function(df, path){
    df <- sapply(df, as.character)
    df[is.na(df)] <- ""
    write.table(as.data.frame(df), path, quote=FALSE, 
            row.names=FALSE, sep="\t", na = "")
}

add_func_annot <- function(df, annots, annot_header = TRUE, annot_idcol = 1){
    if (is.character(annots)){
        annots <- read_tsv(annots, col_names = annot_header, col_types = cols())
    }
    id_name <- names(annots[annot_idcol])    
    full_annot <- left_join(df, annots, by = setNames(id_name, "feature") ) 

    return(full_annot)
}

plot_pk_distr <- function(peak_table, sample){
    data <- peak_table %>% select(width, peak) %>% unique
    p2 <- ggplot(data, aes(x=width)) +
      stat_bin(breaks = seq(0,20000,1000), position="identity") +
      xlim(0, 20000)
    pg <- ggplot_build(p2)
    
    pk_med <- peak_table %>% select(width, peak) %>% unique %>% pull(width) %>% median
    p3 <- ggplot(pg$data[[1]], aes(x=x, y=y, fill = x)) +
      geom_bar(stat="identity", color = "darkgreen", width = 1000) +
      scale_x_continuous("Island size (kb)", expand = c(0.02,0),
                         breaks = seq(0, 20000, by = 3000),
                         labels = seq(0,20, by = 3)) + 
      scale_y_continuous("Number of regions", expand = c(0,0),
                         limits = c(0, 1.1* max(pg$data[[1]]$y)),
                         breaks = seq(0, 1.1* max(pg$data[[1]]$y), by = 500)) +
      theme_classic() +
      theme(legend.position="none",
           plot.margin = margin(t = 5, r = 10, b = 5, l = 7, unit = "pt")) + 
      geom_vline(xintercept = pk_med, color = "blue") + 
      geom_text(aes(x = 8500 + pk_med, label=paste0("Median = ",round(pk_med/1000, 2)," kb"), 
                    y = 4000), colour="blue", size=3.7) +
      scale_fill_gradient(low="lightgreen", high="darkgreen")
    p3
    ggsave(p3, file = paste0("output_figures/peak_len_hist_",sample,".pdf"), width = 2.5, height = 3 )
}

## find results
sample_1 <- args[1]  # treatment
sample_2 <- args[2]  # control
peak_sw <- "epic2"   ## hardcoded

bed_1 <- list.files(".", pattern = paste0(peak_sw,".*",sample_1,".*.bed"), 
                    full.names = TRUE, ignore.case = TRUE)
bed_2 <- list.files(".", pattern = paste0(peak_sw,".*",sample_2,".*.bed"), 
                    full.names = TRUE, ignore.case = TRUE)
manorm <- list.files(pattern = "*MAvalues.xls", recursive = TRUE)


### combine MAnorm and epic2 res
grMN <- read_tsv(manorm, col_types = cols())
grMN <- makeGRangesFromDataFrame(grMN, keep.extra.columns = TRUE)

md_1 <- get_metadat(bed_1, grMN, name = sample_1)
md_2 <- get_metadat(bed_2, grMN, name = sample_2)
values(grMN) <- cbind(md_1, md_2, values(grMN))


## annotate manorm results
gff_file <- "genome_annot.gff"
annot_file <- "genome_annot.txt"
annot_header <- TRUE
annot_idcol <- 1
out_path <- paste0(tools::file_path_sans_ext(manorm), "-annot_ol500.tsv")

txdb <- makeTxDbFromGFF(file=gff_file, format="gff", 
                        dataSource="Bra_v3.0", organism="Brassica rapa")
annoData <- toGRanges(txdb, feature="gene")

peak_table <- annotate_peaks(grMN, annoData)
full_annot <- add_func_annot(peak_table, annot_file, annot_header, annot_idcol) 
write_annot(full_annot, out_path)

### stats, output list of changing genes (filter by P and M)
pval <- 0.05     ## hardcoded
mval <- as.numeric(args[3])

peak_table %>% filter(P_value < pval, M_value >= mval) %>% 
  pull(feature) %>% unique %>% na.omit %>%
  write_lines(paste0(dirname(manorm),"/DMG_P",pval,"-Mge",mval,".txt"))
peak_table %>% filter(P_value < pval, M_value <= -mval) %>% 
  pull(feature) %>% unique %>% na.omit %>%
  write_lines(paste0(dirname(manorm),"/DMG_P",pval,"-Mle",mval,".txt"))

### annotate epic2 results
out_path <- paste0(tools::file_path_sans_ext(bed_1), "-annot_ol500.tsv")
peak_table <- annotate_peaks(bed_1, annoData)
full_annot <- add_func_annot(peak_table, annot_file, annot_header, annot_idcol) 
plot_pk_distr(peak_table, sample_1)
write_annot(full_annot, out_path)

out_path <- paste0(tools::file_path_sans_ext(bed_2), "-annot_ol500.tsv")
peak_table <- annotate_peaks(bed_2, annoData)
full_annot <- add_func_annot(peak_table, annot_file, annot_header, annot_idcol) 
plot_pk_distr(peak_table, sample_2)
write_annot(full_annot, out_path)

