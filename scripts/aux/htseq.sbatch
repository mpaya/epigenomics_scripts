#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:30:00
#SBATCH --job-name=htseq

source config.txt
mkdir -p $1 && cd $1
base=$(basename $1)

sample=$2
idx=1

while [ $idx -eq 1 ] 
do

    ## Run htseq. Count files will have descriptive names.
    htseq-count \
    --format=sam \
    --order=name \
    --stranded=no \
    --type=gene \
    --idattr=ID \
    $sample \
    "$genome_gff" \
    > $(basename ${sample%%.*})${htseq_pt} \
    2> $(basename ${sample%%.*}).htOut
     
    ## If there are errors, prepare files for next rounds. Save original sam as sam.bak
    a=($(grep -l Error *Out))
    
    if [[ -n ${a[@]} ]]
    then
        file=$(grep -l Error *Out)
        suf=bak
        badline=$(grep Error $file | sed 's;.*line \([0-9]*\) .*;\1;' | head -1)
        badfile=$(grep Error $file | sed 's;.*file \(.*\)):;\1;' | head -1)
        badread=$(head -$badline $badfile | tail -1 | cut -f1)
        if [[ -e "$badfile".bak ]]; then suf=bak2; fi
        mv "$badfile" "$badfile"."$suf"
        grep -v $badread $badfile."$suf" > "$badfile" 
        htseq_maxerr=$((htseq_maxerr-1))
        if [ $htseq_maxerr -eq 0]; then echo "$htseq_maxerr errors reached" > errlog; idx=0; fi
    else
        idx=0
    fi
    wait 

done
