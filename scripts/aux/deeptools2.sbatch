#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 5:00:00
#SBATCH --job-name=deept2

source config.txt

cpus=4
sampledir="$1"
origdir="$2"

mkdir -p "$sampledir" && cd "$sampledir"

## get info from metadata
id=${expidx[$(basename ${sampledir%%/analy*})]}
reps=($($(awk -v i=$id -v j=$col_sample \
    'NR>1{if ($i != "") {print $j}}' FS='\t' "$metadata") | sort -u)
muts=($(awk -v i=$id -v j=$col_strain \
    'NR>1{if ($i != "") {print $j}}' FS='\t' "$metadata") | sort -u)

## link files
cnt=0
declare -A filelist
for s in ${reps[@]}
do
    file=$(find "$origdir"/$s -name *$finalbam_pt)
    ln -s "$file" $s.bam
    ln -s "$file".bai $s.bam.bai
    filelist[$cnt]=$s.bam
#    eval arr${muts[$cnt]//[^a-zA-Z0-9]/}\+\=\($s.bam\)
    cnt=$(( cnt + 1 ))
done

echo "files linked"

## multiBamSummary
multiBamSummary bins \
--numberOfProcessors "$cpus"  \
--outFileName 'deeptools_coverage_matrix' \
--bamfiles $(echo ${filelist[@]}) \
--binSize '10000' \
--distanceBetweenBins '0' \
--outRawCounts ChIP_counts.tab \
--scalingFactors ChIP_cnt_sf.tab
# --labels $(echo ${reps[@]})  # this option is bugged

## for each condition with 2 or more samples, create separate matrices
# for s in $(awk '{if (a[$1] == 1) {print $1;}; a[$1]+=1;}' RS=' ' <(echo ${muts[@]}))
# do
#   multiBamSummary bins \
#   --numberOfProcessors "$cpus"  \
#   --outFileName deeptools_cov_$s.matrix \
#   --bamfiles $(eval echo \${arr${s//[^a-zA-Z0-9]/}[@]}) \
#   --binSize '10000' \
#   --distanceBetweenBins '0' \
#   --outRawCounts ChIP_counts_"$s".tab \
#   --scalingFactors ChIP_cnt_sf_"$s".tab
# done

echo "multiBamSummary done"

## plotPCA 
plotPCA \
--corData 'deeptools_coverage_matrix' \
--plotTitle 'ChIP-Seq mutants' \
--plotFile 'deeptools_pca.png' \
--plotFileFormat 'png'

echo "plotPCA done"

## plotCorrelation
plotCorrelation \
--corData 'deeptools_coverage_matrix' \
--plotFile 'deeptools_corr.png' \
--corMethod 'spearman' \
--whatToPlot 'heatmap'  \
--colorMap 'RdYlBu' \
--plotNumbers \
--plotTitle 'ChIP-Seq mutants'  \
--plotWidth 11.0 \
--plotHeight 9.5 \
--skipZeros \
--plotFileFormat 'png'

echo "plotCorrelation done"

## scatterplot
Rscript $scriptd/deeptools2.R 

## cleanup
rm *.bam *.bai
echo "cleanup done"


