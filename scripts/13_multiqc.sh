#!/bin/bash
[[ ! -z $(which Rscript) ]] || module load gcc/6.4.0 R/3.6.0
source config.txt

exp=(rnaseq chipseq)

install=0

if [ $install -eq 1 ]
then
    pip install multiqc==1.7
fi

for f in ${exp[@]}
do
    cd "$basedir"/$f/analysis
    if [ ! -d *multiqc ] 
    then
        step=$(ls -d [0-9]*/ | sed 's;\([0-9]\+\).*;\1;' | sort -un | awk 'END{print $1+1}')
        mkdir "$step"_multiqc
    fi
    cd *multiqc
    ln -s "$metadata" md.txt
    multiqc -f --dirs -dd 2 ../{1..5}*/ &> multiqc.log
    Rscript $scriptd/multiqc.R &>> multiqc.log
    rm md.txt

    ## GATHER EPIC2 RES
	epic_dirs=($(ls -d ../*epic2))
	[ ! -z $epic_dirs ] && 
	echo "## EPIC2 RESULTS ON GAP $epic2_gap WITH FILTER log2FC $peak_fc" > epic2_Res.tsv && \
	for d in ${epic_dirs[@]}
	do
	    echo -e "## Results from directory $d, $(ls -d ../4${d:4:1}* | sed 's;[[:punct:]]\+4.\?_\(.*\);\1;')" >> epic2_Res.tsv
	    sed -n 1p $(ls $d/*/*gapRes.txt | head -1) >> epic2_Res.tsv
	    for s in $d/*; 
	    do 
	        awk -v g=$epic2_gap '{if ($2 == g){print $0}}' $s/*gapRes.txt
	    done | sort -u >> epic2_Res.tsv
	    echo -e "\n\n" >> epic2_Res.tsv
	done

    ## GATHER MANORM RES
	mn_files=($(find ../*manorm -name *all_MAvalues.xls))
	[ ! -z $mn_files ] && \
	echo -e "manorm_dir\tsample1\tsample2\tmapping_sw\ttotal_peaks_1\ttotal_peaks_2\t\
	    unique_peaks_1\tunique_peaks_2\tcommon_merged_peaks\tmanorm_peaks\tm_gt_0\t\
	    m_lt_0\tm_ge_$manorm_m\tm_le_$manorm_m" > manormRes.tsv && \
	for f in ${mn_files[@]}
	do
	    aln=$(ls -d ../4${f:4:1}* | cut -d_ -f2)
	    s1=$(grep "Sample 1 name" "${f%/*}"/manorm_*.log | sed 's;.*= ;;')
	    s2=$(grep "Sample 2 name" "${f%/*}"/manorm_*.log | sed 's;.*= ;;')
	    tpk1=$(grep "Total peaks of sample 1" "${f%/*}"/manorm_*.log | sed 's;.*1: \([0-9]\+\).*;\1;')
	    tpk2=$(grep "Total peaks of sample 2" "${f%/*}"/manorm_*.log | sed 's;.*2: \([0-9]\+\).*;\1;')
	    upk1=$(grep "Total peaks of sample 1" "${f%/*}"/manorm_*.log | sed 's;.*e: \([0-9]\+\).*;\1;')
	    upk2=$(grep "Total peaks of sample 2" "${f%/*}"/manorm_*.log | sed 's;.*e: \([0-9]\+\).*;\1;')
	    cpk=$(grep "merged common" "${f%/*}"/manorm_*.log | sed 's;.*: ;;')
	    mnpk=$(sed 1d "${f}" | wc -l)
	    up=$(awk -F '\t' 'NR>1{m_val=$5; if(m_val>0){a++;}} END{print a}' "${f}")
	    do=$(awk -F '\t' 'NR>1{m_val=$5; if(m_val<0){a++;}} END{print a}' "${f}")
	    filtup=$(awk -F '\t' -v m=$manorm_m 'NR>1{if($5 >= m){a++;}} END{print a}' "${f}")
	    filtdo=$(awk -F '\t' -v m=$manorm_m 'NR>1{if($5 <= -m){a++;}} END{print a}' "${f}")
	    
	    echo -e "${f%/*}\t$s1\t$s2\t$aln\t$tpk1\t$tpk2\t$upk1\t$upk2\t$cpk\t$mnpk\t\
	        $up\t$do\t$filtup\t$filtdo" >> manormRes.tsv
	done

done

