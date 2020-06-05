#!/bin/bash
for f in logs/jobs*.log;
do 
    [ -d ${f%%.log} ] || mkdir -p ${f%%.log}; 
    for id in $(cat $f | grep -v '##' | cut -d' ' -f4); 
    do 
        [ -f slurm-$id.* ] && mv slurm-$id.* ${f%%.log} # 2> /dev/null
    done
done
