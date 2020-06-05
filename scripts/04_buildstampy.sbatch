#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 3:00:00
#SBATCH --job-name=ref_idx

module load python/2.7.15
module load gcc/6.4.0 bwa/0.7.17 samtools/1.9
source config.txt
install=0

if [ $install -eq 1 ]
then
	[ -d "$sw_dir" ] || mkdir -p "$sw_dir" && cd "$sw_dir"
	
	wget https://005.medsci.ox.ac.uk/files-library/stampy-latest.tgz
	tar xzvf stampy-latest.tgz
	ln -s stampy-1.0.32/stampy.py ~/bin/
	
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5/bowtie2-2.3.5-linux-x86_64.zip
	unzip bowtie2-2.3.5-linux-x86_64.zip 
	ln -s bowtie2-2.3.5-linux-x86_64/bowtie2* ~/bin/
    
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
	tar -jxvf bwa-0.7.17.tar.bz2 && cd bwa-0.7.17/ && make
	ln -s bwa ~/bin/
    
    conda install -c bioconda samtools=1.9
    
    cd $sw_dir && wget https://github.com/broadinstitute/picard/releases/download/2.21.1/picard.jar
    ln -s picard.jar ~/bin
fi

[ -d "$refdir" ] || mkdir -p "$refdir" && cd "$refdir"

## Build a genome (.stidx) file:
stampy.py --species="$genome_sp" --assembly="$genome_id" \
            -G "$genome_id" "$genome_file" &> log_stampy-idx.out

## Build a hash (.sthash) file:
stampy.py -g "$genome_id" -H "$genome_id" &> log_stampy-hash.out

## BWA index
bwa index -p "$genome_id" "$genome_file" &> log_bwa-idx.out

## Bowtie2 build
bowtie2-build "$genome_file" "$genome_id" &> log_bt2-build.out
