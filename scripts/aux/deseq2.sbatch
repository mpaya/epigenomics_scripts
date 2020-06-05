#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --job-name=deseq2

source config.txt
module load gcc/6.4.0 R/3.6.0
mkdir -p $1 && cd $1

sampledir=$1
htseqdir=$2
pattern=$3

[ -f anno.dat ] || ln -s "$metadata" anno.dat
[ -d counts/ ] || mkdir counts
[ -f genome_annot.txt ] || ln -s "$genome_annot" genome_annot.txt

files=($(find "$htseqdir" -name $htseq_pt))
for f in ${files[@]}
do
    [ -f counts/${f##*/} ] || ln -s "$f" ./counts/ 
done

## run main DESeq2 script
Rscript $scriptd/deseq2.R

## run DESeq2 on each batch
Rscript $scriptd/deseqByBatch.R

