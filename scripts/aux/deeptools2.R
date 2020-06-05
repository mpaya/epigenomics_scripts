#!/usr/bin/Rscript

if (!requireNamespace("tidyverse", quietly = TRUE)){
    install.packages("tidyverse", repos='http://cran.us.r-project.org')}
if (!requireNamespace("GGally", quietly = TRUE)){
    install.packages("GGally", repos='http://cran.us.r-project.org')}
if (!requireNamespace("fields", quietly = TRUE)){
    install.packages("fields", repos='http://cran.us.r-project.org')}

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

plot_chipres <- function(cnt, cnt_sf, out, m = "") {
    df <- read_tsv(cnt, col_types = cols()) %>% as.data.frame
    rownames(df) <- df %>% select(1:3) %>% unite(bin) %>% pull
    colnames(df) <- df %>% colnames %>% str_replace_all("'","") %>% str_replace_all(".bam","")
    df <- df %>% select(4:ncol(df))
    
    # subset solumns by higher sample
    df <- df %>% select(names(df) %>% str_subset(pattern = m) %>% sort)
    
    df_sf <- read.table(cnt_sf, header = T, row.names = 1)
    
    df <- sapply(colnames(df), function(x) as.matrix(df)[rownames(df), x] * df_sf[x,]) %>%
      as.data.frame(., row.names = rownames(.)) %>%
      filter(rowSums(.) > ncol(.))
    
    p <- GGally::ggpairs(log2(df+1),
                        1:ncol(df),
                        lower = list(continuous = GGscatterPlot),
                        upper = list(continuous = wrap("cor", method= "spearman")))
    ggsave(p, file = out, width = ncol(df) + 1, height = ncol(df) + 1, compression="lzw")
    return(p)
}

### ChIP-Seq

## plot results from all samples
plot_chipres(cnt = "ChIP_counts.tab",cnt_sf = "ChIP_cnt_sf.tab", out = "corr_chipseq.tiff")

## scatterplot replicates from individual matrices
# muts <- list.files(pattern = "ChIP_counts_.*[.]tab") %>% str_replace_all("ChIP_counts_(.*)[.]tab","\\1")
# for (m in muts){
#     plot_chipres(cnt = paste0("ChIP_counts_",m,".tab"),cnt_sf = paste0("ChIP_cnt_sf_",m,".tab"), 
#                  out = paste0("corr_chipseq_",m,".tiff"))
# }

## scatterplot replicates, find strains from dir names and subset full matrix if replicated
muts <- list.dirs(recursive = FALSE) %>% str_remove_all("./")
for (m in muts){
    if (length(list.files(pattern = paste0(m, ".*.bam"))) > 1) {
        plot_chipres(cnt = "ChIP_counts.tab",cnt_sf = "ChIP_cnt_sf.tab", 
            out = paste0("corr_chipseq_",m,".tiff"), m)
    }
}
