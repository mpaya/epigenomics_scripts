#!/usr/bin/Rscript
library(tidyverse)

## data files
md <- read_tsv("md.txt", col_types = cols()) %>% select(full_sample, sample)
f1 <- "multiqc_skewer.txt"
f2 <- "multiqc_picard_AlignmentSummaryMetrics.txt"
f3 <- "multiqc_bowtie2.txt"
f4 <- "multiqc_htseq.txt"
f5 <- "multiqc_picard_dups.txt"
f6 <- "multiqc_fastqc.txt"

## functions
gather_res <- function(df){
  r1 <- df %>% group_by(.dots = c("sample","prog")) %>% select(-starts_with("pct")) %>% summarise_if(is.numeric, sum)
  r2 <- r1 %>% group_by(prog) %>% summarise_if(is.numeric, sum) %>% add_column(`.before` = 1, sample = "total / average")
  
  p1 <- df %>% group_by(.dots = c("sample","prog")) %>% select(sample,prog,starts_with("pct_")) %>% summarise_if(is.numeric, mean)
  p2 <- p1 %>% group_by(prog) %>% summarise_if(is.numeric, mean) %>% add_column(`.before` = 1, sample = "total / average")
  p3 <- p1 %>% group_by(prog) %>% summarise_if(is.numeric, sd) %>% add_column(`.before` = 1, sample = "sd")
  
  stats <- full_join(bind_rows(r1,r2),bind_rows(p1,p2,p3), by = c("sample", "prog"))
  return(stats)
}

## core code
dd <- dir(pattern = "[^.]*data$")
for (d in dd){
  setwd(dirname(d))
  ## raw fastqc
  if (file.exists(paste0(d, "/", f6))) {
    df1 <- read_tsv(list.files(d, f6, full.names = T), col_types = cols()) %>% 
      separate(Sample, c("step","full_sample","filename"), " [|] ") %>% 
      filter(step == (.) %>% pull(step) %>% head(1)) %>%
      mutate(full_sample = str_remove_all(full_sample, "_[12]$"), pct_len = as.numeric(`Sequence length`)) %>% 
      select(full_sample, raw_reads = `Total Sequences`, pct_len, pct_GC = `%GC`) %>% 
      right_join(md, ., by = "full_sample") %>% unique %>% add_column(prog = "read_qual")
  }
  
  ## read preprocessing
  if (file.exists(paste0(d, "/", f1))) {
    df1 <- read_tsv(list.files(d, f1, full.names = T), col_types = cols()) %>% 
      separate(Sample, c("step","full_sample","filename"), " [|] ") %>%
      select(full_sample, total_reads = r_avail, pct_avail, pct_trimmed, pct_untrimmed) %>%
      right_join(df1, ., by = "full_sample")
  }
  ## gather read results
  read_stats <- df1 %>% gather_res() %>% select(-prog) %>% rename(raw_read_len = pct_len)
  df1 %>% rename(raw_read_len = pct_len) %>% write_tsv("1_read_fullRes.tsv")
  write_tsv(read_stats, "1_read_stats.tsv")
  
  ## read mapping
  if (file.exists(paste0(d, "/", f2))) {
    write_full_res <- 1
    # picard on mapped and filtered alignments
      df2 <- read_tsv(list.files(d, f2, full.names = T), col_types = cols()) %>% 
        separate(Sample, c("step","full_sample","file"), " [|] ") %>%
        select(full_sample, step, total_reads = TOTAL_READS, total_mapped = PF_READS_ALIGNED, mapped_hq = PF_HQ_ALIGNED_READS, 
               pct_mapped = PCT_PF_READS_ALIGNED, pct_mismatch = PF_MISMATCH_RATE, read_len = MEAN_READ_LENGTH) %>%
        mutate_at(vars(total_reads, total_mapped, mapped_hq), list(~ ./2)) %>% 
        right_join(md, ., by = "full_sample") %>% 
        separate(step, sep = c(1,2,3,20), into = c("step","sub","s","prog"))
      
      ## Add trimmed read lengths to read report.
      df1 <- df2 %>% dplyr::filter(grepl("^4",step)) %>% select(c(ends_with("sample"), pct_read_len = read_len)) %>% 
        unique%>% left_join(df1, ., by = c("sample","full_sample"))
      read_stats <- df1 %>% gather_res() %>% select(-prog) %>% rename(trim_read_len = pct_read_len, raw_read_len = pct_len)
      df1 %>% rename(trim_read_len = pct_read_len, raw_read_len = pct_len) %>% write_tsv("1_read_fullRes.tsv")
      write_tsv(read_stats, "1_read_stats.tsv")
      
      ## mapping results
      map_res <- df2 %>% dplyr::filter(grepl("^4",step)) %>% 
        mutate(pct_map_hq = mapped_hq / total_reads) %>% 
        mutate(unmapped = total_reads - total_mapped, pct_unmap = unmapped/total_reads) %>% 
        select(c(contains("sample"), sub, prog, total_reads, contains("map"), starts_with("pct")))
      
      # add bowtie2
      if (file.exists(paste0(d, "/", f3))) {
        map_res <- read_tsv(list.files(d, f3, full.names = T), col_types = cols()) %>% 
          separate(Sample, c("step","full_sample","file"), " [|] ") %>% 
          separate(step, "_", into = c("pref","prog")) %>%
          select(full_sample, prog, total_reads, multimap = paired_aligned_multi) %>%
          full_join(map_res, ., by = c("full_sample", "prog", "total_reads")) %>%
          mutate(pct_multi = multimap/total_reads)
      }
      
      ## add ChIP-Seq deduplication
      if (file.exists(paste0(d, "/", f5))) {
        df5 <- read_tsv(list.files(d, f5, full.names = T), col_types = cols()) %>% 
          separate(Sample, c("step","full_sample","filename"), " [|] ") %>%
          separate(step, sep = c(1,2,3,20), into = c("step","sub","s","prog"))
        map_res <- df5 %>% filter(grepl("filter$",prog)) %>% 
          select(full_sample, total_dupes = READ_PAIR_DUPLICATES, sub) %>% 
          full_join(df2 %>% filter(step == "4") %>% select(full_sample, sub, total_reads), ., by = c("full_sample","sub")) %>%
          mutate(pct_dupes = total_dupes/total_reads) %>% 
          full_join(map_res, ., by = c("full_sample","sub","total_reads"))
        
      }
      
      ## add ChIP-Seq filtering steps
      if (nrow(filter(df2, grepl("filter$",prog))) > 0){
        map_res <- df2 %>% filter(grepl("filter$",prog)) %>% 
          select(full_sample, filt_reads = total_reads, sub) %>% 
          full_join(df2 %>% filter(step == "4") %>% select(full_sample, sub, total_reads), ., by = c("full_sample","sub")) %>%
          mutate(pct_filt = filt_reads/total_reads) %>% 
          full_join(map_res, ., by = c("full_sample","sub","total_reads"))
      }
      
      ## add RNA-Seq counts
      if (file.exists(paste0(d, "/", f4))) {
        df4 <- read_tsv(list.files(d, f4, full.names = T), col_types = cols()) %>% 
          separate(Sample, c("step","full_sample","file"), " [|] ") %>%
          separate(step, sep = c(1,2,3,20), into = c("step","sub","s","prog2")) %>% 
          full_join(df2 %>% select(sub, prog, full_sample, sample) %>% unique, ., by = c("sub","full_sample")) %>% 
          select(-file, -s, -prog2, -sub, -step, total_reads = total_count, pct_assigned = percent_assigned)
        
        ## gather count results
        count_stats <- gather_res(df4)
        write_tsv(df4, "3_count_fullRes.tsv")
        write_tsv(count_stats, "3_count_stats.tsv")
        
        df4b <- df4 %>% select(-pct_assigned) %>% rename_if(is.numeric, list(~ paste0("pct_",.)))
        cat("\nNormalization with mapped reads\n", file = "3_count_stats.tsv", append = TRUE)
        df4b %>% mutate(pct_mapped = pct_total_reads - pct_not_aligned) %>% 
          mutate_if(is.numeric, list(~ ./pct_mapped*100)) %>% gather_res() %>% 
          write_tsv("3_count_stats.tsv", append = TRUE, col_names = TRUE)
        
        map_res <- df4 %>%
          mutate(pct_assigned = pct_assigned/100) %>% 
          select(full_sample, prog, cnt_assigned = assigned, pct_assigned) %>% 
          full_join(map_res, ., by = c("full_sample","prog"))
      }
        
      ## gather final alignment results
      map_res <- map_res %>% mutate_at(vars(starts_with("pct")), list(~ .*100))
      map_stats <- gather_res(map_res)
      write_tsv(map_res, "2_map_fullRes.tsv")
      write_tsv(map_stats, "2_map_stats.tsv")
  }
}

if (write_full_res == 1){
  map_res %>%
    full_join(df1 %>% select(-prog,-total_reads), . , by = c("full_sample", "sample")) %>% 
    mutate_at(vars(starts_with("pct")), list(~ ./100)) %>%
    rename(Sample = full_sample) %>% 
    mutate_at(vars(starts_with("pct")), scales::percent) %>% 
    mutate("Total PE reads" = raw_reads, 
           "Reads processed" = paste0(total_reads,"\n(",pct_avail,")"),
           "Reads mapped" = paste0(total_mapped,"\n(",pct_mapped,")")) %>% 
    mutate("Reads duplicated" = if (exists('total_dupes')) paste0(total_dupes,"\n(",pct_dupes,")") else NA) %>% 
    mutate("Multi-mapped reads" = if (exists('multimap')) paste0(multimap,"\n(",pct_multi,")") else NA) %>% 
    mutate("Non-mapping reads" = paste0(unmapped,"\n(",pct_unmap,")")) %>% 
    arrange(prog, Sample) %>%
    select(prog, matches("^[A-Z]", ignore.case = FALSE)) %>%
    write_tsv("Mapping_Res.tsv")
}

