#!/bin/bash
# https://devhints.io/bash Bash scripting cheatsheet -> replace first match (1 ->2 ) remove prefix, suffix ...
# Tons of stackoverflow like this one https://stackoverflow.com/questions/47302898/redirect-stdout-and-stderr-to-file-permanently-but-keep-printing-them
# script to get the results from final_project

# for conda activate to work
# https://stackoverflow.com/questions/59665231/conda-4-7-how-to-activate-env-in-bash-script
# 1st Create a log file with all stdout and stderr from the script keeping 
# sourcing bash_profile (needed for conda initialization)
# If sourced from main script initialization is done. No need to source .bash_profile
source ~/scripts/init_gdav $1


# 2. Genome Analysis (basic checks)
srcpath=~/final_project/RNAseq
dstpath=~/results/genome_analysis/basic_checks
mkdir -p $dstpath

echo -e "\n"
echo "      ########     genome analysis (basic checks)     ######## "
echo -e "\n"

echo check directory and files
ls -la $srcpath


# we have 8 files, we asume 4 samples (each one with 2 files, forward and reverse), 2 for high temperature and 2 for normal temperature 
# First check if the files are fastq format, 4 lines for entry. check for 5 entries (4*5 = 20 lines)
echo -e "\n"
echo cat first 20 lines to check file format 
echo
zcat $srcpath/hightemp01.r1.fq.gz | head -n 20
echo
echo fastq format, 4 lines each entry. It seems all sequences have maximum quality in all elements \"I\"
echo -e "\n"


echo "Check for samples (r1 and r2)"
echo

nsamples=0
for f_file in $srcpath/*.r1.fq.gz
do 
    r_file=${f_file/r1/r2}
    #motus profile -f metagenomics-hotspring-hightemp.1.fq.gz -r metagenomics-hotspring-hightemp.2.fq.gz -o ../results/metagenomics/hightemp.motu -n hightemp
    # remove prefix
    sample=${f_file##$srcpath/}
    # remove suffix
    sample=${sample%%.r1.fq.gz}
    echo "first line from sample $sample (r1 and r2)"
    #echo zcat $f_file | head -n 1
    zcat $f_file | head -n 1
    #echo
    #echo zcat $r_file | head -n 1
    zcat $r_file | head -n 1
    ((nsamples+=1))
    echo
done

echo "We have $nsamples samples"
echo

# All sequence reads have the name ^@AQUIFEX
echo "checking numbers of reads (counts line 1 of each entry \"^@AQUIFEX_\" ) and max quality (all \"I\")"
echo
for f_file in $srcpath/*.r1.fq.gz
do 
    r_file=${f_file/r1/r2}
    #motus profile -f metagenomics-hotspring-hightemp.1.fq.gz -r metagenomics-hotspring-hightemp.2.fq.gz -o ../results/metagenomics/hightemp.motu -n hightemp
    # remove prefix
    sample=${f_file##$srcpath/}
    # remove suffix
    sample=${sample%%.r1.fq.gz}
    echo -n "Reads from sample $sample (r1 / r2) and all \"I\" quality :  "
    #echo zcat $f_file | head -n 1
    f_reads=$(zcat $f_file | grep "^@AQUIFEX_" | wc -l)
    r_reads=$(zcat $r_file | grep "^@AQUIFEX_" | wc -l)

    f_maxq_reads=$(zcat $f_file | grep -E "^I[I]+I$" | wc -l)
    r_maxq_reads=$(zcat $r_file | grep -E "^I[I]+I$" | wc -l)

    echo "r1: $f_reads / r2: $r_reads max \"I\" quality r1: $f_maxq_reads / r2: $r_maxq_reads"
    if [ $f_reads == $r_reads ] && [ $f_reads == $f_maxq_reads ] && [ $f_maxq_reads == $r_maxq_reads ] ; then
        echo "  - forward and reverse have same reads and all elements in all reads have \"I\" (40) quality (99.99% probability of the element to be correct)"
    fi
done

echo
echo -n "All sequences have the same length: "
zcat $srcpath/hightemp01.r1.fq.gz | head -n 4 | grep -E "^I[I]+I$" | wc -L


