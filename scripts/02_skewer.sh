#!/bin/bash

## Script to run skewer on CESGA
## It requires a set of	fastq files in provided	dirs and a metadata file
## first metadata column is sample names,         
## second column is plant line
## other metadata collumns include prefix ($idx) and endings ($idx + 1) of read pairs
## original fastq are linked to	homogenize sample names	across experiments (skewer.sh)

source config.txt
script=aux/skewer.sh
stepdir=2_skewer
install=0

# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_02.log

if [ $install -eq 1 ]
then
	## download file with sequencing adapters
	[ -d "$sw_dir" ] || mkdir -p "$sw_dir" && cd "$sw_dir"
	wget https://gist.githubusercontent.com/photocyte/3edd9401d0b13476e60f8b104c2575f8/raw/7e1e7e3dac674c196d8a05169558b675a57b7e23/Sequencing_adaptors.fasta
	## Get skewer:
	[ -d ~/bin/ ] || mkdir ~/bin/
	[ -d "$sw_dir" ] || mkdir -p "$sw_dir" && cd "$sw_dir"
	wget https://sourceforge.net/projects/skewer/files/Binaries/skewer-0.2.2-linux-x86_64
	chmod a+x skewer-0.2.2-linux-x86_64
	ln -s skewer-0.2.2-linux-x86_64 ~/bin/skewer
fi

for e in ${exp[@]}
do
	readdir="$basedir"/$e/ARCHIVE
	sample_rows=($(awk -v i=${expidx[$e]} 'NR>1{if ($i != "") {print NR}}' FS='\t' "$metadata"))
	for idx in ${sample_rows[@]}
	do
	    files=($(ls "$readdir"/$(sed -n "$idx"p "$metadata" | cut -f ${expidx[$e]}-$((${expidx[$e]} + 1)) | tr '\t' '*')))
		samplename=$(sed -n "$idx"p "$metadata" | cut -f $col_sample)
		sampledir="$basedir"/$e/analysis/$stepdir/"$samplename"
		$command $script "$sampledir" ${files[0]} ${files[1]}  >> logs/jobs_02.log
	done
done

