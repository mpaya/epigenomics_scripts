#!/bin/bash

## Script to run fastqc	on CESGA
## It requires a set of	PE fastq files in provided	dirs and a metadata file
## first metadata column is sample names,         
## second column is plant line
## other metadata columns include prefix ($idx) and endings ($idx + 1) of read pairs
## original fastq are linked to	homogenize sample names	across experiments (fastqc.sh)

source config.txt
script=aux/fastqc.sbatch
stepdir=1_fastqc
install=0

# add option to run locally
command=$1
if [ $command == "" ]; then command=sbatch; fi

# code to get message from prompt
message="new batch submission"
echo "# $message"  >> logs/jobs_01.log

## install anaconda3 if needed
# cd $sw_dir && \
#     wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh && \
#     bash Anaconda3-2019.10-Linux-x86_64.sh -b -p $anaconda_dir

if [ $install -eq 1 ]
then
    # install fastqc on a custom python environment
    conda create --name fastqc python=3.7 fastqc=0.11.8 -c bioconda
fi

for e in ${exp[@]}
do
	readdir="$basedir"/$e/ARCHIVE
	sample_rows=($(awk -v i=${expidx[$e]} 'NR>1{if ($i != "") {print NR}}' FS='\t' "$metadata"))
	for idx in ${sample_rows[@]}
	do
	    files=($(for i in $(sed -n "$idx"p "$metadata" | cut -f ${expidx[$e]}-$((${expidx[$e]} + 1))); 
                 do ls "$readdir"/$i; done))
		samplename=$(sed -n "$idx"p "$metadata" | cut -f $col_sample)
		sampledir="$basedir"/$e/analysis/$stepdir/"$samplename"
		# a suffix is added to include compression info (e.g. .fq.gz)
	    suf=."${files[0]#*.}"
		$command $script "$sampledir" ${files[0]} ${files[1]} $suf >> logs/jobs_01.log
	done
done

