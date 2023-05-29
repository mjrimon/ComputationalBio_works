#!/bin/bash
# https://devhints.io/bash Bash scripting cheatsheet -> replace first match (1 ->2 ) remove prefix, suffix ...
# Tons of stackoverflow like this one https://stackoverflow.com/questions/47302898/redirect-stdout-and-stderr-to-file-permanently-but-keep-printing-them
# script to get the results from final_project

# for conda activate to work
# https://stackoverflow.com/questions/59665231/conda-4-7-how-to-activate-env-in-bash-script
# we have to run this script in interactive mode 'bash -i' or source ~/.bash_profile

# 1st Create a log file with all stdout and stderr from the script keeping 
# sourcing bash_profile (needed for conda initialization)
# If sourced from main script initialization is done. No need to source .bash_profile

source ~/scripts/init_gdav $1


# 4. Genome Analysis (variant calling)

gensrcpath=~/final_project
srcpath=~/results/genome_analysis/read_mapping
dstpath=~/results/genome_analysis/variant_calling
mkdir -p $dstpath

echo index reference genome
genomefile=$srcpath/aquifex_genome.fasta


echo "Activating bcftools environment"
echo

conda activate bcftools

echo "show stats for merged file"
bcftools stats $dstpath/merged.vcf | grep -P "ST\t|SN\t|IDD\t" | grep -vP ":\t0"
echo

echo "show variants with quality > 100"
bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/' | awk -F '\t' '$6 > 100 {print}'
echo
echo -n "Total variants with quality > 100: "
bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^#" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/' | awk -F '\t' '$6 > 100 {print}' | wc -l
echo 

echo "show variants with depth of coverage > 100"
bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk -F '\t' '$8 > 100  {print}'
echo
echo -n "Total variants with depth of coverage > 100: "
bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^#" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk -F '\t' '$8 > 100  {print}' | wc -l
echo


echo "Comparing the output from merged BAM file whith the separate files"
echo 
echo "merged BAM file"
bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk -F '\t' '{print $2 "\t" $4 "\t" $5 "\t" $6 "\t" $8}'
echo

echo "separate BAM files"
bcftools view $dstpath/variant_calling_separately.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk -F '\t' '{print $2 "\t" $4 "\t" $5 "\t" $6 "\t" $8}'
echo


echo "4.5 Identify the variant with best quality"
echo

genomegff=$gensrcpath/genome.gff


export LC_NUMERIC="en_US.UTF-8" 

bestvariant=$(bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^#" | tr ";" "\t" | awk -F '\t' ' {print $2 "\t" $4 "\t" $5 "\t" $6 "\t" $8}' \
                                                                  |  sort -t$'\t' -k4 -nr | awk -F '\t' 'NR==1 {print $1}')


cat $genomegff | awk -v ref="$bestvariant" -F '\t' '{ if ((ref <= $5 ) && (ref >= $4 )){print}}'
echo 




echo 4.6.- 


echo "Activating samtools environment"
echo

conda activate samtools

echo Merge by temperature sorted bam files into one and index for use in IGV
samtools merge -h $dstpath/h*1.sorted.bam $dstpath/hightemp_merged.bam  $dstpath/h*.sorted.bam

samtools merge -h $dstpath/n*1.sorted.bam $dstpath/normaltemp_merged.bam  $dstpath/n*.sorted.bam

samtools index $dstpath/hightemp_merged.bam
samtools index $dstpath/normaltemp_merged.bam

#bcftools view variant_calling_separately.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk -F '\t' '{print $2 $4 $5 $6 $8}'


#bcftools view $dstpath/variant_calling_merged.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'| awk '$6 > 50 && $8 > 10  {print}' | wc -l
#bcftools view variant_calling_merged.bcf | grep -v "^##" | tr ";" "\t" | sed 's/DP=//g'| awk '$6 > 50 && $8 > 10  {print}' | wc -l
#| grep -v "^##" | tr ";" "\t" | sed 's/DP=//g;s/INFO/DP/'



#bcftools view variant_calling_merged.bcf | grep -v "^#" | awk -F '\t' '$6 > 100 {print}' | wc -l

