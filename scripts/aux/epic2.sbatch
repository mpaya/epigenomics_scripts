#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --job-name=epic2
#SBATCH --mem-per-cpu=15GB

module load gcc/6.4.0 bedtools/2.28.0 R/3.6.0
source "$(pwd)"/config.txt  # "${SLURM_SUBMIT_DIR}"
cpus=1

samples=$2
wt_f=$3
input_f=$4
 
compare_wt=1
compare_reps=0

## sample vs wt
if [ $compare_wt -eq 1 ] && [[ "$wt_name" != ${1##*/} ]]
then
	mkdir -p "$1"_vs_"$wt_name" && cd "$1"_vs_"$wt_name"
    # if analysis failed and needs re-run, move files to cwd
    [[ -d "tracks" ]] && mv */* .

	$scriptd/1a_epic2_gaps_analysis.sh \
	-c $2 \
	-i $input_f \
	-C $3 \
	-I $input_f \
	-1 ${1##*/}  \
	-2 $wt_name \
	-f $genome_file \
	-e $egf \
	-g $genome_gff \
	-M 4 \
	-L $filterlist \
	-t $cpus \
	-F $peak_fc
	
	#cd ..
fi

## compare duplicates

if [ $compare_reps -eq 1 ]
then
	mkdir -p "$1" && cd "$1"
	a=($(echo $samples | tr ',' '\n'))
	[ ${#a[@]} -eq 2 ] && $scriptd/1a_epic2_gaps_analysis.sh \
	-c ${a[0]} \
	-i $4 \
	-C ${a[1]} \
	-I $4 \
	-1 $(basename ${a[0]%%.*}) \
	-2 $(basename ${a[1]%%.*}) \
	-f $genome_file \
	-e $egf \
	-g $genome_gff \
	-M 4 \
	-L $filterlist \
	-r ${1##*/} \
	-t $cpus
	
	cd ..
fi

