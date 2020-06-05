#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5:00:00
#SBATCH --job-name=deept1

cpus=1

source config.txt
module load samtools/1.9 picard/2.18.14
PATH=$EBROOTPICARD/:$PATH
mkdir -p $1 && cd $1

## link bam files
ln -s $2 .
ln -s "$2".bai .

## if two files, link and merge them
if [ ! -z $3 ]; 
then
    ln -s $3 .
    ln -s "$3".bai .

    ## MergeSamFiles
    java -Xmx2g -jar picard.jar \
    MergeSamFiles  \
    INPUT="$2"  \
    INPUT="$3"  \
    OUTPUT="${1##*/}".merged.bam \
    MERGE_SEQUENCE_DICTIONARIES="false"  \
    ASSUME_SORTED="false"  \
    USE_THREADING=true \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY="LENIENT" 

    samtools index "${1##*/}".merged.bam
fi

## calc egs
genome_len=$(echo $(( $(grep -v '>' $genome_file | wc -c) - $(grep -v '>' $genome_file | wc -l) )))
if [ ! -z $egf ]
then
    egs=$(echo $genome_len $egf | awk '{printf "%4.0f",$1*$2}')
else
    egs=$genome_len
fi

## bamCoverage
for f in *.bam
do
    bamCoverage \
    --numberOfProcessors "${cpus}"  \
    --bam "$f" \
    --outFileName $(sed 's/[a-z]\+.bam/bigwig/' <(echo "${f}")) \
    --outFileFormat 'bigwig'  \
    --binSize 50  \
    --normalizeUsing RPGC \
    --effectiveGenomeSize $egs
done

