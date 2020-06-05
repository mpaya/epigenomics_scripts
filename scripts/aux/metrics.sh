#!/bin/bash

source ${0%/aux*}/config.txt
base=$1
cpus=$2

## convert to indexed sorted bam
samtools view \
-b $base.sam \
-h \
-@ $(( cpus - 1 )) \
-o $base.bam

samtools sort \
-o $base${sortbam_pt} \
-@ $(( cpus - 1 )) \
$base.bam

rm -f $base.bam
samtools index $base${sortbam_pt}

## calculate metrics
java -Xmx2g -jar $EBROOTPICARD/picard.jar \
CollectAlignmentSummaryMetrics \
VALIDATION_STRINGENCY=LENIENT \
R="$genome_file" \
I=$base${sortbam_pt} \
O=$base.metrics 2> picard.err

