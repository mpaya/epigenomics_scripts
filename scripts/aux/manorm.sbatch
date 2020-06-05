#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5:00:00
#SBATCH --job-name=manorm
#SBATCH --mem-per-cpu=10GB

module load gcc/6.4.0 samtools/1.9 picard/2.18.14 R/3.6.0
PATH=$EBROOTPICARD/:$PATH
source config.txt

outdir="$1"
bamdir="$2"
beddir="$3"
s=${bamdir##*/}
pattern_bed=*"$s"_"$epic2_gap".bed

cpus=1
compare_batch=1

mkdir -p $1 && cd $1

## find files of sample 1
# if merged file was created, first find it
# else, get original files
# merge of bam files is done in tmp folder to avoid clash among jobs

bamfiles1=($([[ -d tmp_"$s" ]] && find tmp_"$s" -name "$s"*merged.bam))
if [ ${#bamfiles1[@]} -eq 0 ];
then
    bamfiles1=($(find "$bamdir"* -name *$dup_pt))
fi
bedfile1=($(find "$beddir"* -name $pattern_bed))

## if two files, merge them
if [ ${#bamfiles1[@]} -eq 1 ]; 
then
    bamfile1=${bamfiles1[0]}
elif [ ${#bamfiles1[@]} -eq 2 ]
then
   echo merging $s
    [[ -d tmp_"$s" ]] || mkdir tmp_"$s"
    bambase=tmp_"$s"/"$s".dupmark
    ## MergeSamFiles
    java -Xmx8g -jar picard.jar \
    MergeSamFiles  \
    INPUT="${bamfiles1[0]}"  \
    INPUT="${bamfiles1[1]}"  \
    OUTPUT="$bambase".merged.bam \
    MERGE_SEQUENCE_DICTIONARIES="false"  \
    ASSUME_SORTED="false"  \
    USE_THREADING=true \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY="LENIENT" 2> $bambase.mergelog
    samtools index "$bambase".merged.bam
    bamfile1="$bambase".merged.bam
elif [[ ${#bamfiles2[@]} -gt 2 ]]
then
    echo "More than two replicates on $s not supported" && k=1
fi
echo "sample1 files processed"

## find files for wild type
wt_beddir=${beddir/$s/*$wt_name}

if [[ compare_batch -eq 1 ]]
then
   f=$(basename $(pwd | sed 's;analy.*;;'))
   batch=$(awk -v i=${expidx[$f]} -v j=$col_strain -v b=$col_batch -v s=$s \
      'NR>1{if ($i != "" && $j == s) {print $b}}' FS='\t' "$metadata" | sort -u)
   ## do something if length batch > 1
   wt_name=$(awk -v i=${expidx[$f]} -v j=$col_strain -v k=$col_sample \
      -v wt=$wt_name -v b=$col_batch -v sb=$batch 'NR>1{if ($i != "" && 
      $j == wt && $b == sb) {print $k}}' FS='\t' "$metadata")
   ## same, what if multiple...
fi

wt_bamdir=${bamdir/$s/$wt_name}

bamfiles2=($([[ -d tmp_"$s" ]] && find tmp_"$s" -name "$wt_name"*merged.bam))
if [ ${#bamfiles2[@]} -eq 0 ];
then
    bamfiles2=($(find "$wt_bamdir"* -name *$dup_pt))
fi
bedfile2=($(find "$wt_beddir"* -name ${pattern_bed/$s/$wt_name}))

## if two files, merge them
if [ ${#bamfiles2[@]} -eq 1 ]; 
then
    bamfile2=${bamfiles2[0]}
elif [ ${#bamfiles2[@]} -eq 2 ]
then
   echo merging $wt_name
    [[ -d tmp_"$s" ]] || mkdir tmp_"$s"
    bambase=tmp_"$s"/"$wt_name".dupmark
    ## MergeSamFiles
    java -Xmx2g -jar picard.jar \
    MergeSamFiles  \
    INPUT="${bamfiles2[0]}"  \
    INPUT="${bamfiles2[1]}"  \
    OUTPUT="$bambase".merged.bam \
    MERGE_SEQUENCE_DICTIONARIES="false"  \
    ASSUME_SORTED="false"  \
    USE_THREADING=true \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY="LENIENT" 2> $bambase.mergelog
    samtools index "$bambase".merged.bam
    bamfile2="$bambase".merged.bam
elif [[ ${#bamfiles2[@]} -gt 2 ]]
then
    echo "More than two replicates on $wt_name not supported" && k=1
fi
echo "sample2 files processed"

if [[ ${#bamfiles1[@]} -eq 0 ]]; then echo "missing bam files for $s"; k=1; fi
if [[ ${#bedfile1[@]} -eq 0 ]]; then echo "missing bed file for $s"; k=1; fi
if [[ ${#bamfiles2[@]} -eq 0 ]]; then echo "missing bam files for $wt_name"; k=1; fi
if [[ ${#bedfile2[@]} -eq 0 ]]; then echo "missing bed file for $wt_name"; k=1; fi
if [[ $k -eq 1 ]]; then exit 1; fi

## run MAnorm
manorm_dir=manorm_"$s"_vs_"$wt_name"

## run manorm if no previous run exists
if [ ! -f manorm_"$s"_vs_"$wt_name"/*xls ] 
then
    source $anaconda_dir/etc/profile.d/conda.sh
    conda activate $anaconda_dir/envs/manorm
    manorm \
    --p1 ${bedfile1} \
    --p2 ${bedfile2} \
    --peak-format bed \
    --r1 $bamfile1 \
    --r2 $bamfile2 \
    --read-format bam \
    --name1 $s \
    --name2 $wt_name \
    --paired-end \
    -o "$manorm_dir"/ \
    2> manorm_"$s".log
    conda deactivate
    ## add numbers to log
    echo -e "\n# peaks\tM>0\tM>0.1\tM>0.25\tM>0.5\tM>1" >> manorm_"$s".log
    awk -F '\t' 'NR>1{m_val=sqrt($5^2); if(m_val>0){a++;} if(m_val>.1){b++;} 
        if(m_val>.25){c++;} if(m_val>.5){d++;} if(m_val>1){e++;} }END{print "total",a,b,c,d,e}' \
        OFS='\t' "$manorm_dir"/*xls >> manorm_"$s".log
    awk -F '\t' 'NR>1&&$5>0{m_val=$5; if(m_val>0){a++;} if(m_val>.1){b++;} 
        if(m_val>.25){c++;}if(m_val>.5){d++;} if(m_val>1){e++;} }END{print "M > 0",a,b,c,d,e}' \
        OFS='\t' "$manorm_dir"/*xls >> manorm_"$s".log
    awk -F '\t' 'NR>1&&$5<0{m_val=-$5; if(m_val>0){a++;} if(m_val>.1){b++;} 
        if(m_val>.25){c++;} if(m_val>.5){d++;} if(m_val>1){e++;} }END{print "M < 0",a,b,c,d,e}' \
        OFS='\t' "$manorm_dir"/*xls >> manorm_"$s".log
    mv manorm_"$s".log manorm_"$s"_vs_"$wt_name"
    ## finish
    echo "manorm done"
else
    echo "using previous manorm run"
fi

## run annotation
cd "$manorm_dir"
ln -s ${bedfile1} .
ln -s ${bedfile2} .
ln -s "$genome_gff" genome_annot.gff
ln -s "$genome_annot" genome_annot.txt

Rscript $scriptd/manorm.R $s "$wt_name" $manorm_m

rm *bed genome_annot*

echo "manorm analysis finished"
