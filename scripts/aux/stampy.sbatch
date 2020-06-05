#!/bin/bash
#SBATCH -n 1
#SBATCH -c 3
#SBATCH -t 18:00:00
#SBATCH --job-name=stampy
#SBATCH --mem-per-cpu=10GB

module load python/2.7.15
module load gcc/6.4.0 bwa/0.7.17 samtools/1.9 picard/2.18.14
source config.txt
cpus=8

mkdir -p $1 && cd $1
base=$(basename $1)
ref="$refdir"/$genome_id

## To speed up mapping, use BWA (recommended) and multithreading:

bwa aln -q10 -t"$cpus" $ref $2 > "$base"_1.sai 2>> bwa.log
bwa aln -q10 -t"$cpus" $ref $3 > "$base"_2.sai 2>> bwa.log
bwa sampe $ref "$base"_1.sai "$base"_2.sai $2 $3 2>> bwa.log | \
    samtools view -Sb - > bwa.bam 2>> bwa.log

stampy.py \
-g $ref \
-h $ref \
-t "$cpus" \
--readgroup=ID:$base,SM:$base,LB:lib1,PL:illumina \
-o $base.sam \
--bamkeepgoodreads \
-M bwa.bam 2> $base.err

bash "$scriptd"/metrics.sh $base $cpus

