#!/bin/bash

## Script to run fastqc after skewer on CESGA
## It requires a set of	fastq files in provided	dirs and a metadata file
## first metadata column is sample names,         
## second column is plant line
## other metadata collumns include prefix ($idx) and endings ($idx + 1) of read pairs
## original fastq are linked to	homogenize sample names	across experiments (fastqc2.sh)

source config.txt
script=aux/fastqc.sh
stepdir=3_fastqc
sourcedir=2_skewer
skewer_pt=-pair?.fastq

# add option to run locally
command=$1
if [ $command == "" ]; then command=sbatch; fi

# code to get message from prompt
message="new batch submission"
echo "# $message"  >> logs/jobs_03.log

for i in $(seq 0 $(( ${#exp[@]} - 1 )))
do
	readdir="$basedir"/${exp[$i]}/analysis/$sourcedir
	sample_rows=($(awk -v i=${expidx[${exp[$i]}]} 'NR>1{if ($i != "") {print NR}}' FS='\t' "$metadata"))
	for idx in ${sample_rows[@]}
	do
    	files=($(ls "$readdir"/$(awk -v i=$idx 'NR==i{print $1}' "$metadata")/*$skewer_pt))
		samplename=$(sed -n "$idx"p "$metadata" | cut -f $col_sample)
		sampledir="$basedir"/${exp[$i]}/analysis/$stepdir/"$samplename"
    	suf=-trimmed."${files[0]#*.}"
		$command $script "$sampledir" ${files[0]} ${files[1]} $suf  >> logs/jobs_03.log
	done
done
