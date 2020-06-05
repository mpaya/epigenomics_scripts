#!/bin/bash
source config.txt
script1=aux/epic2.sbatch
script2=aux/bins.sbatch
stepdir=(6b_epic2 6a_epic2)
sources=(5b_filter 5a_filter)
exp=(chipseq)

## installation
if [[ $install -eq 1 ]]
then
    # $anaconda_dir/bin/conda install -c bioconda epic2
    pip install epic2==0.0.34
    # install bedtools
    cd $sw_dir
    wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
    chmod a+x bedtools
    ln -s bedtools ~/bin
fi

# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_10.log

## submit jobs
for i in $(seq 0 $(( ${#exp[@]} - 1 )))
do
    samples=($(awk -v i=${expidx[${exp[$i]}]} -v j=$col_strain \
        'NR>1{if ($i != "") {print $j}}' FS='\t' "$metadata" | sort -u))
    for x in $(seq 0 $(( ${#sources[@]} - 1 )))
    do
        origdir="$basedir"/${exp[$i]}/analysis/${sources[$x]}
        # stepdir and sources are checked for same length
        if [ ${#stepdir[@]} -eq ${#sources[@]} ]; then y=$x; else y=0; fi
        for s in ${samples[@]}
        do
            files=$(ls "$origdir"/$s*/*$finalbam_pt | paste -sd,)
            input=$(ls "$origdir"/${input_name}*/*$finalbam_pt | paste -sd,)
            wt_f=$(ls "$origdir"/${wt_name}*/*$finalbam_pt | paste -sd,)
            ## run epic2 gap optimization
            sampledir="$basedir"/${exp[$i]}/analysis/${stepdir[$y]}/$s
            $command $script1 "$sampledir" "${files}" "${wt_f}" "${input}" >> logs/jobs_10.log
            ## run epic2 bin optimization
            binsdir="$basedir"/${exp[$i]}/analysis/${stepdir[$y]}_bin/$s
            $command $script2 "$binsdir" "${files}" "${wt_f}" "${input}" >> logs/jobs_10.log
        done
    done
done

## an option to see all samples together is to plot score by gap
## PSEUDOCODE
if [[ $collect -eq 1 ]]
then
    cd "$basedir"/${exp[$i]}/analysis/${stepdir[$i]}/
    Rscript -e 'flow_res <- list();
    res <- list.files(path, pattern = ".*/epic_rawRes/*txt");
    for (idx in seq(length(res))){
        flow_res[[idx]] <- read_tsv(paste(path, res[[idx]], sep = "/"), col_types = cols()) %>% 
            mutate(peak_len = End - Start);
        df <- add_row(df, tissue = "flowers", gap = idx - 1, n_islands = length(flow_res[[idx]]$Score), 
                      mean_len = mean(flow_res[[idx]]$peak_len), score = sum(flow_res[[idx]]$Score));
    }
    df <- df[-1,];
    ggplot() + 
    geom_point(data=df, aes(x = mean_len, y = score, color = tissue, shape = gap), size=5) + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_shape_identity() + 
    labs(x= "Mean length", y = "Aggregated score");
    
    ggsave("epic_gapSize.jpg", width = 5, height = 4);'
fi