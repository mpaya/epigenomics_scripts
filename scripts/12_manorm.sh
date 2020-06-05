#!/bin/bash
source config.txt
script=aux/manorm.sbatch
stepdir=(8a_manorm 8b_manorm)
bamsources=(5a_filter 5b_filter)
bedsources=(6a_epic2 6b_epic2)
exp=(chipseq)

install=0

# add option to run locally
command=$1
if [[ $command == "" ]]; then command=sbatch; fi

# code to get message from prompt
message="Batch submission: $(date --iso-8601=seconds)"
echo "# $message" >> logs/jobs_12.log

if [[ $install -eq 1 ]]
then

    ## Build MAnorm 1.2.0 from source
    cd $sw_dir
    git clone https://github.com/shao-lab/MAnorm.git
    cd MAnorm/

    ## write configuration file
    cat <<EOM > environment.yml
name: manorm
channels:
 - defaults
 - conda-forge
 - bioconda
dependencies:
 - python=3.6.6
 - openssl=1.0.2p=h14c3975_1002
 - pysam=0.10.0
prefix: $anaconda_dir/envs/manorm
EOM

    ## create environment and install
    conda env create -f environment.yml
    conda activate manorm
    pip install .
    conda deactivate
fi


## submit jobs
for f in ${exp[@]}
do
    for i in $(seq 0 $(( ${#stepdir[@]} - 1 )))
    do
        muts=($(awk -v i=${expidx[$f]} -v j=$col_strain \
            'NR>1{if ($i != "") {print $j}}' FS='\t' "$metadata" \
            | grep -v $input_name | grep -v $wt_name | sort -u))
        for s in ${muts[@]};
        do
            bamdir="$basedir"/$f/analysis/${bamsources[$i]}/$s
            beddir="$basedir"/$f/analysis/${bedsources[$i]}/$s
            outdir="$basedir"/$f/analysis/${stepdir[$i]}
            $command $script "$outdir" "$bamdir" "$beddir" >> logs/jobs_12.log
        done
    done
done

