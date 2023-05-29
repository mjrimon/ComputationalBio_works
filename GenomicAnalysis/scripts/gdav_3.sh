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

# 3. Genome Analysis (read mapping)

srcpath=~/final_project/RNAseq
dstpath=~/results/genome_analysis/read_mapping
mkdir -p $dstpath

echo -e "\n"
echo "      ########     genome analysis (read mapping)     ######## "
echo -e "\n"

# File with genome
genomefile=~/final_project/genome.fasta
cp $genomefile $dstpath/aquifex_genome.fasta

# Creating index file from genome
bwa index -p $dstpath/aquifex_genome $genomefile
echo

for f_file in $srcpath/*.r1.fq.gz
do 
    r_file=${f_file/r1/r2}
    #bwa mem aquifex_genome hightemp01.r1.fq.gz hightemp01.r2.fq.gz > hightemp01.sam 2> hightemp01.sam.err

    # remove prefix
    sample=${f_file##$srcpath/}
    # remove suffix
    sample=${sample%%.r1.fq.gz}
    
    #Creating sam file from indexed genome and fastq files
    bwa mem $dstpath/aquifex_genome $f_file $r_file > $dstpath/$sample.sam 2> $dstpath/$sample.sam.err

done
echo -e "\n"

echo "Activating samtools environment"
echo

conda activate samtools

# convert sam files to bam files (binary format) 
for samfile in $dstpath/*.sam
do
    # remove suffix
    outfile=${samfile%%.sam}
    echo Sample ${outfile##$dstpath/}, creating binary bam file from sam file
    #echo samtools view -S -b -h  $samfile -o $outfile.bam
    samtools view -S -b -h  $samfile -o $outfile.bam
    rm $samfile
    echo
    echo Flag stat data from ${outfile##$dstpath/}
    samtools flagstat $outfile.bam | grep read
    echo
done
echo -e "\n"

