#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --job-name=epicbin
#SBATCH --mem-per-cpu=10GB

module load gcc/6.4.0 bedtools/2.28.0 R/3.6.0
source $(pwd)/config.txt
cpus=1

samples=$2
wt_f=$3
input_f=$4
 
compare_wt=1
compare_reps=0

# if analysis failed and needs re-run, move files to cwd
[[ -d "tracks" ]] && mv */* .

if [ $compare_wt -eq 1 ] && [[ "$wt_name" != ${1##*/} ]]
then
	mkdir -p "$1"_vs_"$wt_name" && cd "$1"_vs_"$wt_name"
    $scriptd/1c_epic2_bins_analysis.sh \
    -c $2 \
    -i $input_f \
    -C $3 \
    -I $input_f \
    -1 ${1##*/}  \
    -2 $wt_name \
    -f $genome_file \
    -e $egf \
    -g $genome_gff \
    -L $filterlist \
    -t $cpus
fi
