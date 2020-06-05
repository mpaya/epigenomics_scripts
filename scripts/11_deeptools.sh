#!/bin/bash
source config.txt
stepdir=(7b_deeptools 7a_deeptools)
sources=(5b_filter 5a_filter)
exp=(chipseq)

# stepdir and sources are checked for same length

if [ install -eq 1 ]
then
    if [[ -z $($anaconda_dir/bin/conda config --show | grep bioconda) ]]
    then
        $anaconda_dir/bin/conda config --add channels bioconda
    fi
    $anaconda_dir/bin/conda create --name deeptools python=3.7 deeptools -y
fi
$anaconda_dir/bin/conda activate deeptools

# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_11.log

## merge files and get bigwig
script=aux/deeptools.sbatch
## submit jobs
for f in ${exp[@]}
do
    for idx in $(seq 0 $(( ${#sources[@]} - 1 )))
    do
        sourcedir=${sources[$idx]}
        if [ ${#stepdir[@]} -eq ${#sources[@]} ]; then i=$idx; else i=0; fi
        awk 'NR>1{if ($5 != ""){print $1 FS $2}}' FS='\t' "$basedir"/mutant_list.csv > tmp
        for s in $(cut -f2 tmp | sort -u);
        do
            files=($(ls "$basedir"/$f/analysis/$sourcedir/$s*/$s*$finalbam_pt))
            sampledir="$basedir"/$f/analysis/${stepdir[$i]}/$s
            sbatch $script "$sampledir" ${files[0]} ${files[1]} >> logs/jobs_11.log
        done
    done
done
rm tmp

## plot pca and heatmap
script=aux/deeptools2.sbatch
## submit jobs
for f in ${exp[@]}
do
    for idx in $(seq 0 $(( ${#sources[@]} - 1 )))
    do
        if [ ${#stepdir[@]} -eq ${#sources[@]} ]; then i=$idx; else i=0; fi
        origdir="$basedir"/$f/analysis/${sources[$idx]}
        sampledir="$basedir"/$f/analysis/${stepdir[$i]}
        $command $script "$sampledir" "$origdir" >> logs/jobs_11.log
    done
done

$anaconda_dir/bin/conda deactivate
