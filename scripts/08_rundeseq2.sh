#!/bin/bash
source config.txt
script=aux/deseq2.sbatch
stepdir=(6a_deseq2 6b_deseq2)
sources=(5a_htseq 5b_htseq)
exp=(rnaseq)


# add option to run locally
command=$1
if [ $command == "" ]; then command=sbatch; fi

# code to get message from prompt
message="new batch submission"
echo "# $message" >> logs/jobs_08.log

## submit jobs
for i in $(seq 0 $(( ${#exp[@]} - 1 )))
do
    for x in $(seq 0 $(( ${#sources[@]} - 1 )))
    do
        # stepdir and sources are checked for same length
        if [ ${#stepdir[@]} -eq ${#sources[@]} ]; then y=$x; else y=0; fi
        htseqdir="$basedir"/${exp[$i]}/analysis/${sources[$x]}
        sampledir="$basedir"/${exp[$i]}/analysis/${stepdir[$y]}
        $command $script "$sampledir" "$htseqdir" >> logs/jobs_08.log
    done
done
