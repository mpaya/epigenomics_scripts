#!/bin/bash
#===============================================================================
#
#          FILE:  1a_epic2_gaps_analysis.sh
#
#         USAGE:  ./1a_epic2_gaps_analysis.sh
#
#   DESCRIPTION:
#
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:  Miriam Paya Milans
#       COMPANY:  CBGP
#       VERSION:  1.0
#       CREATED:  06/06/2019
#      REVISION:  ---
#===============================================================================

## Some help for usage
usage_short()
{
cat <<EOM
    Usage: script.sh [-c CHIP1] [-i INPUT1] [-C CHIP2] [-I INPUT2]
                     [-1 NAME1] [-2 NAME2] [-m MIN_GAP] [-M MAX_GAP]
                     [-f FASTA_FILE] [-e EGF] [-g GFF_FILE] [-x] [-F FC]
                     [-L LIST_FILE] [-r REPS_NAME] [-t THREADS] [-h]

EOM
}

usage_long()
{
echo -e '\tThis script runs epic2 over multiple gap sizes and compares results.'
echo -e '\tSoftware requirements:'
echo -e '\t\texport PATH="$epic2":"$bedtools2";"$Rscript":$PATH '
cat <<EOM
Options:
  -m, --min_gap=MIN_GAP             minimum gap size (0)
  -M, --max_gap=MAX_GAP             maximum gap size (7)
  -c, --chip1=CHIP1[,CHIP]          ChIP files for sample 1 (c1_chip*.[bed,bam])
  -i, --input1=INPUT1[,INPUT]       INPUT files for sample 1 (c1_input*.[bed,bam])
  -1, --name1=NAME1                 Name for first sample (col0)
  -C, --chip2=CHIP2[,CHIP]          ChIP files for sample 2 (c2_chip*.[bed,bam])
  -I, --input2=INPUT2[,INPUT]       INPUT files for sample 2 (c1_input*.[bed,bam])
  -2, --name2=NAME2                 Name for second sample (elf6)
  -f, --fasta=FASTA_FILE            FASTA file with genome sequences
  -e, --egf=EGF                     Effective (i.e. mappable) genome fraction (0.99)
  -g, --gff=GFF_FILE                GFF file with gene annotations (*[gff,gff3])
  -L, --chr_list=LIST_FILE          List of chromosomes to keep on coverage calculations
                                      (e.g. discard scaffolds or organelles)
  -F, --peak_fc                     Minimum log2FC on output peaks (0)
  -r, --reps                        Common name for merging replicated samples
  -t, --threads=THREADS             Number of threads to be used (6)
  -x, --force                       Delete all current files before running
  -h                                Display usage and exit
  --help                            Display this help and exit

EOM
}


## set default values
min_gap=0
max_gap=7

name1="sample1"
name2="sample2"
chip_files_1=$(ls c1_chip* 2>/dev/null)
input_files_1=$(ls c1_input* 2>/dev/null)
chip_files_2=$(ls c2_chip* 2>/dev/null)
input_files_2=$(ls c2_input* 2>/dev/null)

egf=0.99
cpus=6
reps=""

# get data from arguments
while [ "$1" != "" ]; do
    case $1 in
        -m | --min_gap )        shift
                                min_gap=$1
                                ;;
        -M | --max_gap )        shift
                                max_gap=$1
                                ;;
        -1 | --name1 )          shift
                                name1=$1
                                ;;
        -2 | --name2 )          shift
                                name2=$1
                                ;;
        -c | --chip1 )          shift
                                chip_files_1=($(sed 's;,; ;g' <(echo $1)))
                                ;;
        -i | --input1 )         shift
                                input_files_1=($(sed 's;,; ;g' <(echo $1)))
                                ;;
        -C | --chip2 )          shift
                                chip_files_2=($(sed 's;,; ;g' <(echo $1)))
                                ;;
        -I | --input2 )         shift
                                input_files_2=($(sed 's;,; ;g' <(echo $1)))
                                ;;
        -f | --fasta )          shift
                                genome_fasta=$1
                                ;;
        -e | --egf )            shift
                                egf=$1
                                ;;
        -g | --gff )            shift
                                gff=$1
                                ;;
        -L | --chr_list )       shift
                                fl=$1
                                ;;
        -F | --peak_fc )        shift
                                min_fc=$1
                                ;;
        -r | --reps )           shift
                                reps=$1
                                ;;
        -t | --threads )        shift
                                cpus=$1
                                ;;
        -x | --force )          shift
                                force=TRUE
                                ;;
        -h )                    usage_short
                                exit
                                ;;
        --help )                usage_short
                                usage_long
                                exit
    esac
    shift
done

namelist=($name1 $name2 $reps)

# if FORCE is activated, delete old files
[ ! -z $force ] && echo "delete all" && rm -fr *

## check requirements
# in input data
if [ $min_gap -gt $max_gap ]
then
   echo "min gap must be smaller than max gap"
   exit 1
fi

if [[ ! -z $genome_fasta ]]; 
then
    ln -sf $genome_fasta genome.fa
elif [[ ! -f genome.fa ]]
then
    genome_fasta=($(ls *fa *fasta 2>/dev/null))
    [[ ${#genome_fasta[@]} -eq 1 ]] && ln -sf $genome_fasta genome.fa
    if [[ -z $genome_fasta || ${#genome_fasta[@]} -ne 1 ]];
    then
        echo "unable to determine genome file"
        exit 1
    fi
fi

if [[ ! -z $gff ]]; 
then
    ln -sf $gff annot.gff
elif [[ ! -f annot.gff ]]
then
    gff=($(ls *gff *gff3 2>/dev/null))
    [[ ${#gff[@]} -eq 1 ]] && ln -sf $gff annot.gff
    if [[ -z $gff || ${#gff[@]} -ne 1 ]];
    then
        echo "unable to determine GFF file"
        exit 1
    fi
fi

# in PATH
epic_path=$(which epic2)
bedtools_path=$(which bedtools)
if [[ -z $bedtools_path || -z $epic_path ]]
then 
    echo 'bedtools or epic2 not found in $PATH'
    exit 1
fi

if [[ -z $(which Rscript) ]]
then 
    echo 'Rscript not found in $PATH'
    exit 1
fi

## prepare list of chromosomes to use on coverage
[ -z $fl ] && $(grep '>' $genome_file | sed 's;>;;') > chrlist || ln -s $fl chrlist


#################
### Run epic2 ###
#################
awk '/^>/ {if (seqlen) print seqlen;printf substr($1,2) OFS;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' OFS='\t' ${genome_fasta} > chromsizes
n_gaps=$(( $max_gap - $min_gap + 1 ))

seq $min_gap $max_gap | xargs -n 1 -P $cpus -I {} sh -c "epic2 \
-t ${chip_files_1} \
-c ${input_files_1} \
--chromsizes chromsizes \
--effective-genome-fraction ${egf} \
--gaps-allowed {} \
--keep-duplicates \
> epic2_res_${name1}_{}.txt 2> epic2_${name1}_g{}.log" 

if [[ $(ls epic2_res_${name1}*txt | wc -l) -ne $n_gaps ]]; then 
	echo "epic2 on ${name1} failed"
	exit 1
fi

echo "epic2 on files 1 complete"

seq $min_gap $max_gap | xargs -n 1 -P $cpus -I {} sh -c "epic2 \
-t ${chip_files_2} \
-c ${input_files_2} \
--chromsizes chromsizes \
--effective-genome-fraction ${egf} \
--gaps-allowed {} \
--keep-duplicates \
> epic2_res_${name2}_{}.txt 2> epic2_${name2}_g{}.log" 

if [[ $(ls epic2_res_${name2}*txt | wc -l) -ne $n_gaps ]]; then 
	echo "epic2 on ${name2} failed"
	exit 1
fi

echo "epic2 on files 2 complete"

if [[ $reps != "" ]]
then
    seq $min_gap $max_gap | xargs -n 1 -P $cpus -I {} sh -c "epic2 \
        -t ${chip_files_1} ${chip_files_2} \
        -c ${input_files_2} \
        --chromsizes chromsizes \
        --effective-genome-fraction ${egf} \
        --gaps-allowed {} \
        > epic2_res_${reps}_{}.txt 2> epic2_${reps}_g{}.log"
    echo "epic2 on merged reps complete" 
fi


## convert to BED6+3
for file in epic2_res*txt
do 
    awk -v fc=$min_fc 'NR>1{if ($10>=fc){if ($4==0) {pv=500;qv=500}else{pv=-log($4)/log(10); \
      qv=-log($9)/log(10)}; print $1,$2,$3,"island_"NR-1,int($5),$6,$10,pv,qv}}' \
      OFS="\t" ${file} > ${file/txt/bed}
done
wait
echo "conversion to BED complete"

##########################
### jaccard similarity ###
##########################

echo -e "$name1" "\t" "$name2" "\tintersection\tunion-intersection\tjaccard\tn_intersections" > gap_jaccard.txt
for idx in $(seq $min_gap $max_gap)
do
    for i in $(seq $min_gap $max_gap)
    do
        bedtools jaccard \
        -a <(sort -k1,1 -k2n epic2_res_"$name1"_"$idx".bed | cut -f1-4) \
        -b <(sort -k1,1 -k2n epic2_res_"$name2"_"$i".bed | cut -f1-4) | \
        sed 1d | awk -v idx=$idx -v i=$i '{print idx,i,$0}' OFS='\t' >> gap_jaccard.txt
    done
done

if [[ $(cat gap_jaccard.txt | wc -l) -ne $(( $n_gaps * $n_gaps + 1 )) ]]; then 
	echo "jaccard analysis failed"
	exit 1
fi

echo "jaccard similarity complete"

################
### coverage ###
################
cov_counter=0

# gff
echo '## coverage over all genes' > gap_coverage.txt && \
cov_counter=$(( $cov_counter + 1 ))
for n in ${namelist[@]}
do
    for idx in $(seq $min_gap $max_gap)
    do
        cut -f1-4 epic2_res_"$n"_"$idx".bed | coverageBed \
          -a <(awk '$3=="gene"' $gff) -b stdin | \
          awk -v n=$n -v idx=$idx 'END { print "cov.genes",n,idx,s,g,s/g } \
          { s += $11; g += $12; }' FS='\t' OFS='\t' >> gap_coverage.txt && \
          cov_counter=$(( $cov_counter + 1 ))
    done
done

# genome
echo '## genome coverage' >> gap_coverage.txt && \
cov_counter=$(( $cov_counter + 1 ))

## bed of filtered genome
awk 'NR==FNR{a[$1];next;}
  /^>/ {if (seqlen != 0) print seqlen; i=0;
  if (substr($1,2) in a) {printf substr($1,2) "\t1\t";i=1;} seqlen=0; next} \
  i==1{seqlen+=length($0)}END{if (i==1) print seqlen}' \
  chrlist "${genome_fasta}" > genome_bed
  
# ## bed of entire genome
# awk '/^>/ {if (seqlen) print seqlen;printf substr($1,2) "\t1\t";seqlen=0;next} \
#         {seqlen+=length($0)}END{print seqlen}' "${genome_fasta}" > genome_bed

for n in ${namelist[@]}
do
    for idx in $(seq $min_gap $max_gap)
    do
        cut -f1-4 epic2_res_"$n"_"$idx".bed | coverageBed \
          -a genome_bed -b stdin | \
          awk -v n=$n -v idx=$idx 'END { print "cov.genome",n,idx,s,g,s/g } \
          { s += $5; g += $6; }' FS='\t' OFS='\t' >> gap_coverage.txt && \
          cov_counter=$(( $cov_counter + 1 ))
    done
done

# TEs
if [[ $(awk '$3 ~ /transposable_element/' $gff | wc -l ) -gt 0 ]]
then
    echo '## coverage TEs' >> gap_coverage.txt && \
    cov_counter=$(( $cov_counter + 1 ))
    for n in ${namelist[@]}
    do
        for idx in $(seq $min_gap $max_gap)
        do
            cut -f1-4 epic2_res_"$n"_"$idx".bed | coverageBed \
              -a <(awk '$3 ~ /transposable_element/' $gff) -b stdin | \
              awk -v n=$n -v idx=$idx 'END { print "cov.TEs",n,idx,s,g,s/g } \
              { s += $11; g += $12; }' FS='\t' OFS='\t' >> gap_coverage.txt && \
              cov_counter=$(( $cov_counter + 1 ))
        done
    done
fi

if [[ $(cat gap_coverage.txt | wc -l) -ne $cov_counter ]]; then 
	echo "coverage analysis failed"
	exit 1
fi

echo "coverage analysis complete"

############################
### annotation and plots ###
############################

Rscript "$(dirname "$0")/1b_gap_analysis.R"

echo "annotation and visualization complete"

###############
### cleanup ###
###############
mkdir logs tracks epic_rawRes 2>/dev/null
mv *.log logs/
mv *.bed tracks/
mv epic2*.txt epic_rawRes

rm chromsizes genome_bed chrlist
rm genome.fa annot.gff 2>/dev/null

echo "epic2 gap analysis finished"

