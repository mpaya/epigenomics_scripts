#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH --job-name=fastqc

source config.txt
module load fastqc/0.11.7
# conda activate $anaconda3_dir/envs/fastqc

outdir="$1"
suf=$4
mkdir -p "$outdir"_{1,2}

R1="$outdir"_1/$(basename "$outdir")_1"$suf"
R2="$outdir"_2/$(basename "$outdir")_2"$suf"
ln -s "$2" "$R1"
ln -s "$3" "$R2"

for i in 1 2;
do
    fastqc -o "$1"_$i -f fastq $(eval echo \$R$i) > "$1"_$i/log.out 2> "$1"_$i/log.err
done

# conda deactivate