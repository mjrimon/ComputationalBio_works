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


# 5. Differential expression analysis

gensrcpath=~/final_project
genomegff=$gensrcpath/genome.gff

srcpath=~/results/genome_analysis/variant_calling
dstpath=~/results/genome_analysis/DEA
mkdir -p $dstpath

# we have to ensure always the same order. Used to build counts matrix
sampleslist=("hightemp01 hightemp02 normal01 normal02")
samples=($sampleslist)

outfileext=count


echo 5. Differential expression analysis 


echo "Activating samtools environment"
echo

conda activate samtools

for file in $srcpath/*.sorted.bam
do 
    # remove suffix
    outfile=${file%%.bam}
    #Creating index file from bam file
    echo creating index for $file ...
    samtools index $file > $outfile.bam.bai
done
echo -e "\n"


echo "Activating base environment"
echo

conda activate base

for file in $srcpath/*.sorted.bam
do 
    # remove prefix
    outfile=${file##$srcpath/}
    # remove suffix
    outfile=${outfile%%.sorted.bam}
    #counting file
    echo count file: $file ...
    htseq-count -i locus_tag -t CDS -c $dstpath/$outfile.$outfileext $file $genomegff 
done
echo -e "\n"


# we need to load R anyway, we can join the files in bash with the next lines and use the function DESeqDataSetFromMatrix from package DESeq2 
# or we can use the function DESeqDataSetFromHTSeqCount from de same package and use directly the files just created from htseq-count 

# We have to remember we need to ensure same order in the join of files and conditions definition. It's easier to ensure this order if the joining
# of the files and the condition definions are in the same file. Therefore all this logic will be in the R file. Script ended, no more proccess todo

#exit 0 


# keep this note and the next lines to base this decision 
#join  $dstpath/${samples[0]}.$outfileext $dstpath/${samples[1]}.$outfileext | join - $dstpath/${samples[2]}.$outfileext | join - $dstpath/${samples[3]}.$outfileext > $dstpath/raw_samples.${outfileext}s

# minimal filtering sumcounts > 10
#grep -v "^__" $dstpath/raw_samples.${outfileext}s | awk '{if (($2+$3+$4+$5>10)){print}}' > $dstpath/filtered_samples.${outfileext}s
#echo  "Gene $sampleslist" | cat - $dstpath/filtered_samples.${outfileext}s > $dstpath/samples.${outfileext}s

# tab separated columns
#sed -i 's/ /\t/g' $dstpath/samples.counts 


