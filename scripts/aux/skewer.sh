#!/bin/bash
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 10:00:00
#SBATCH --job-name=skewer
#SBATCH --mem-per-cpu=3GB

source config.txt

# create out dir and link fastq files
outdir=$1
mkdir -p "$outdir"

R1="$outdir"/$(basename "$outdir")_1."${2#*.}"
R2="$outdir"/$(basename "$outdir")_2."${3#*.}"
ln -s "$2" "$R1"
ln -s "$3" "$R2"

# run skewer, installed system-wide
skewer \
-x "$adapters_file" \
-m pe \
-q 15 \
-l 30 \
-o "$outdir"/$(basename "$outdir") \
-t 2 \
"$R1" \
"$R2" \
>& "$outdir"/$(basename "$outdir").trimlog

