#!/bin/bash
module load samtools/1.9
source config.txt
script=aux/bowtie2.sbatch
stepdir=4a_bowtie2
sourcedir=2_skewer
skewer_pt=-pair?.fastq

# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_06.log

## submit jobs
for i in $(seq 0 $(( ${#exp[@]} - 1 )))
do
	readdir="$basedir"/${exp[$i]}/analysis/$sourcedir
	sample_rows=($(awk -v i=${expidx[${exp[$i]}]} 'NR>1{if ($i != "") {print NR}}' FS='\t' "$metadata"))
	for idx in ${sample_rows[@]}
	do
		samplename=$(sed -n "$idx"p "$metadata" | cut -f $col_sample)
    	files=($(ls "$readdir"/$samplename/*$skewer_pt))
		sampledir="$basedir"/${exp[$i]}/analysis/$stepdir/"$samplename"
		$command $script "$sampledir" ${files[0]} ${files[1]} >> logs/jobs_06.log
	done
done

# grep -L 'overall alignment' "$basedir"/*seq/analysis/*bowtie2/*/*.stats
