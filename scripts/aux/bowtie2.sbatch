#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 7:00:00
#SBATCH --job-name=bowtie2
#SBATCH --mem-per-cpu=10GB

module load samtools/1.9 picard/2.18.14
source config.txt
cpus=4

mkdir -p $1 && cd $1
base=$(basename $1)
ref="$refdir"/$genome_id

bowtie2 \
-p $cpus \
-x $ref \
-1 $2 \
-2 $3 \
-N 1 \
--no-mixed \
--no-discordant \
--rg-id "$base" \
--rg "SM:$base" \
--rg "PL:ILLUMINA" \
--rg "LB:$base" \
--rg "PU:01" \
--very-sensitive \
-S $base.sam \
2> "$base".stats 

bash "$scriptd"/metrics.sh $base $cpus
