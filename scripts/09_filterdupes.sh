#!/bin/bash
source config.txt
script=aux/filter.sbatch
stepdir=(5a_filter 5b_filter)
sources=(4a_bowtie2 4b_stampy)
exp=(chipseq)


# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_09.log

## submit jobs
for i in $(seq 0 $(( ${#exp[@]} - 1 )))
do
    sample_rows=($(awk -v i=${expidx[${exp[$i]}]} 'NR>1{if ($i != "") {print NR}}' FS='\t' "$metadata"))
    for x in $(seq 0 $(( ${#sources[@]} - 1 )))
    do
        origdir="$basedir"/${exp[$i]}/analysis/${sources[$x]}
        # stepdir and sources are checked for same length
        if [ ${#stepdir[@]} -eq ${#sources[@]} ]; then y=$x; else y=0; fi
        for idx in ${sample_rows[@]}
        do
            samplename=$(sed -n "$idx"p "$metadata" | cut -f $col_sample)
            files=($(ls "$origdir"/$samplename/*$sortbam_pt))
            sampledir="$basedir"/${exp[$i]}/analysis/$stepdir[$y]/"$samplename"
            $command $script "$sampledir" ${files[0]} >> logs/jobs_09.log
        done
    done
done
