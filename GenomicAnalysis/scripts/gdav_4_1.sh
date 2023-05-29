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


# 4. Genome Analysis (variant calling)

gensrcpath=~/final_project
srcpath=~/results/genome_analysis/read_mapping
dstpath=~/results/genome_analysis/variant_calling
mkdir -p $dstpath

echo "Activating samtools environment"
echo

conda activate samtools
for file in $srcpath/*.bam
do
    # remove prefix
    outfile=${file##$srcpath/}
    # remove suffix
    outfile=${outfile%%.bam}
    echo sorting bam file $outfile
    #echo samtools sort -o  $dstpath/$outfile.sorted.bam $file
    samtools sort -o  $dstpath/$outfile.sorted.bam $file
done

echo -e "\n"


echo Merge all sorted bam files into one. We use the header from the last bam file sorted
samtools merge -f -h $dstpath/$outfile.sorted.bam $dstpath/merged.bam  $dstpath/*.sorted.bam

echo index reference genome
genomefile=$srcpath/aquifex_genome.fasta

samtools faidx $genomefile

echo "Activating bcftools environment"
echo

conda activate bcftools

echo  Perform variant calling from merged bam

bcftools mpileup -f $genomefile -o $dstpath/merged.vcf $dstpath/merged.bam
bcftools call -mv -Ob -o $dstpath/variant_calling_merged.bcf $dstpath/merged.vcf

# We repeat everything separately
echo  Perform variant calling from merged bam

bcftools mpileup -f $genomefile -o $dstpath/separately.vcf $dstpath/*.sorted.bam 
bcftools call -mv -Ob -o $dstpath/variant_calling_separately.bcf $dstpath/separately.vcf
echo

