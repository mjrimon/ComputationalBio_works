#!/bin/bash
# script to get the results from final_project

# 1st Create a log file with all stdout and stderr from the script keeping 
# sourcing bash_profile (needed for conda initialization)
# If sourced from main script initialization is done. No need to source .bash_profile

source ~/scripts/init_gdav $1


# 6. Functional analysis

proteomepath=~/final_project
proteomefile=$proteomepath/proteome.faa

srcpath=~/results/genome_analysis/DEA
dstpath=~/results/genome_analysis/functional_analysis
mkdir -p $dstpath

overexpgenesfile=$srcpath/Annotated_Condition_hightemp_vs_normal_result_padj_0.01.csv

echo -e "\n"
echo 6. Functional analysis
echo -e "\n\n"

echo "Extracting sequences of each over expressed gene out from $overexpgenesfile file"
echo -e "\n"

conda activate base

i=0

while read -r line || [[ -n "${line}" ]]; do
    ((i++))
    # remove header file
    [[ "$i" == '1' ]] && continue
    # skip lines with pAdj <> 0
    [[ $(echo $line | cut -d';' -f5 | tr -d '"') != '0' ]] && continue
    genid=$(echo $line | cut -d';' -f2 | tr -d '"')
    sed -n '/'$genid'/,/^>/p' $proteomefile | sed '${/^>/d}' > $dstpath/$genid.faa
    echo extracted sequence from gene $genid
done < $overexpgenesfile

echo -e "\n\n"


