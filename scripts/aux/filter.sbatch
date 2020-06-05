#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --job-name=filter
#SBATCH --mem-per-cpu=20GB

module load samtools/1.9 picard/2.18.14
PATH=$EBROOTPICARD/:$PATH
source config.txt

mkdir -p $1 && cd $1
base=$(basename "$1")
sample=$2

# ## a file with chr to keep
# if [ ! -f "$filterlist" ]
# then
# 	grep '>' "$genome_file" | sed 's;>;;' > "$filterlist"
# fi

## filter/sort mapping results
[[ ! -f $base.filter.bam || ! -f $base.final.bam ]] && \
samtools view \
-F 1804 -q 1 -b "$sample" \
$([ ! -z "$filterlist" ] && cat "$filterlist" | tr '\n' ' ') | \
samtools sort > $base.filter.bam
samtools index $base.filter.bam; 

## mark dupes but do not remove
java -Xmx20g -jar picard.jar \
MarkDuplicates \
I=$base.filter.bam O=$base${dup_pt} M=$base.dup.qc \
VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
2> picardD.err

## remove dupes
samtools view \
-F 1804 -b $base${dup_pt} > $base${finalbam_pt};
samtools index $base${finalbam_pt}

## collect metrics of final bam
java -Xmx20g -jar $picard.jar \
CollectAlignmentSummaryMetrics \
VALIDATION_STRINGENCY=LENIENT \
R="$genome_file" \
I=$base${finalbam_pt} O=$base${finalbam_pt/bam/metrics} \
2> picardM.err

