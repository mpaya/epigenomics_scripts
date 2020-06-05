#!/bin/bash
source config.txt
script=aux/htseq.sbatch
stepdir=(5a_htseq 5b_htseq)
sources=(4a_stampy 4b_bowtie2)
exp=(rnaseq)
install=0

if [ install -eq 1 ]
then
    module load zlib/1.2.11 gcccore/6.4.0 bzip2/1.0.6
    module load gcccore/6.4.0 curl/7.61.1 openssl/1.1.0i
    $anaconda_dir/bin/pip install htseq==0.11.2 
fi

# add option to run locally
command=$1
if [ $command == "" ]; then command=sbatch; fi

# code to get message from prompt
message="new batch submission"
echo "# $message"  >> logs/jobs_07.log

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
            files=($(ls "$origdir"/$samplename/*.sam))
            sampledir="$basedir"/${exp[$i]}/analysis/$stepdir[$y]/"$samplename"
            $command $script "$sampledir" ${files[0]}  >> logs/jobs_07.log
        done
    done
done
