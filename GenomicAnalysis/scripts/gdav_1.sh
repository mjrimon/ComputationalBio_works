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


# 1.- metagenomics

srcpath=~/final_project
dstpath=~/results/metagenomics
mkdir -p $dstpath

timeout=10
continuemsg="Press return or wait for timeout to continue..."

# activating conda environment
conda activate base

# Files names
# metagenomics-hotspring-hightemp.1.fq.gz
echo -e "Process sample files\n"

for f_file in $srcpath/metagenomics-hotspring-*.1.fq.gz
do 
    r_file=${f_file/1/2}
    #motus profile -f metagenomics-hotspring-hightemp.1.fq.gz -r metagenomics-hotspring-hightemp.2.fq.gz -o ../results/metagenomics/hightemp.motu -n hightemp
    # remove prefix
    sample=${f_file##$srcpath/metagenomics-hotspring-}
    # remove suffix
    sample=${sample%%.1.fq.gz}
    
    motus profile -f $f_file -r $r_file -o $dstpath/$sample.motu -n $sample
    if [ $sample == 'hightemp' ]
    then
        for gens in 5 8 10 
        do
            motus profile -f $f_file -r $r_file -o $dstpath/${sample}_g$gens.motu -n ${sample}_g$gens -g $gens
        done 
    fi
done

# Las abundancias relativas están expresadas en numeros decimales.
# El formato de los numeros decimales esta en ingles, utilizan el punto como separador decimal
# para poder operar matematica y logicamente con estos valores necesitamos configurar los numeros
# con formato ingles 
export LC_NUMERIC="en_US.UTF-8"

echo -e "\n"
echo "      ########     Metagenomics     ######## "
echo -e "\n"

#1.1.- and 1.1.1.- 
# The most abundant organism in high-temperature
echo 1.1.- and 1.1.1.- 
echo The most abundant organism in high-temperature
echo

# Presentacion 1
echo -n "The most abundant organism in high-temperature is: "
tail -n +4 $dstpath/hightemp.motu | sort -t$'\t' -k2 -nr | head -n 1
echo -e "\n"

# Presentacion 2
ht_most_ab_org_name=$(tail -n+4 $dstpath/hightemp.motu | sort -t$'\t' -k2 -nr | head -n 1 | cut -f 1)
ht_most_ab_org_percent=$(tail -n+4 $dstpath/hightemp.motu | sort -t$'\t' -k2 -nr | head -n 1 |  awk -F '\t' '{print $2*100"%"}')
echo "Pretty printing with %" 
echo -n "The most abundant organism in high-temperature is: $ht_most_ab_org_name"
echo ", with a relative abundance of:  $ht_most_ab_org_percent"
echo -e "\n"


#1.1.6.- Do you still detect the same organism if you run mOTUs with a more stringent threshold 
# for the number of marker-gene detections (at least 5 genes)? 
# What is the estimated relative abundance in such a case?
echo 1.1.6.- Do you still detect the same organism if you run mOTUs with a more stringent threshold? 
#1.1.6.- Se puede generar el comando aquí para un solo numero de genes o realizarlo durante el loop inicial,
# Ahora esta en los 2 sitios, hay que borrar uno de ellos...


# obteniendo en el inicio todos los archivos necesarios para contestar las preguntas
# aqui generamos el comando para -g 5

    f_file=$srcpath/metagenomics-hotspring-hightemp.1.fq.gz
    r_file=${f_file/1/2}
    #motus profile -f metagenomics-hotspring-hightemp.1.fq.gz -r metagenomics-hotspring-hightemp.2.fq.gz -o ../results/metagenomics/hightemp.motu -n hightemp
    # remove prefix
    sample=${f_file##$srcpath/metagenomics-hotspring-}
    # remove suffix
    sample=${sample%%.1.fq.gz}
    echo motus profile -f $f_file -r $r_file -o $dstpath/${sample}_g5.motu -n ${sample}_g5 -g 5

echo -e "\n"

# si tenemos 3 valores de gens 5 8 10
echo "Alpha diversity in high temp with marker-gene detections for several marker genes"
for gens in 5 8 10 
do
    echo -n "Alpha diversity in high temp with marker-gene detections of $gens are: "
    awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/hightemp_g$gens.motu | wc -l
    echo -n "The most abundant organism in high-temperature with marker-gene detections of $gens is: "
    sort -t$'\t' -k2 -nr $dstpath/hightemp_g$gens.motu | head -n 1
    echo
done
echo -e "\n"

#1.2.- the most abundant organisms in normal-temperature > 1%
# the last field must be greater than 1%  ($NF>0.01) and remove not defined ($1!="-1")
echo "1.2 The most abundant organism in normal-temperature (>1%) are: "
tail -n +4 $dstpath/normaltemp.motu | sort -t$'\t' -k2 -nr | awk '{if ($NF>0.01 && $1!="-1") {print}}'
echo

#1.2.1.- Are they also present in the high-temperature samples? At which relative abundance?

echo "1.2.1.- They are also present in high-temperature sample with relative abundance of:"

most_normal=$(tail -n +4 $dstpath/normaltemp.motu | sort -t$'\t' -k2 -nr | awk -F '\t' '{if ($NF>0.01 && $1!="-1") {print $1}}' | awk '{print $NF}' | sed 's/\[//;s/\]//')

grep "$most_normal" $dstpath/hightemp.motu | sort -t$'\t' -k2 -nr

echo -e "\n"

#1.4.- In your opinion, which is the condition (high or normal temperature) with a greater level of alpha biodiversity?

# alpha diversity, number of organisms detected
# the last field must be greater than 0 ($NF>0) and remove not defined ($1!="-1")
echo "1.4.- which is the condition (high or normal temperature) with a greater level of alpha biodiversity?"
echo -n "Alpha diversity in high temp is: "
awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/hightemp.motu | wc -l

echo -n "Alpha diversity in normal temp is: "
awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/normaltemp.motu | wc -l

echo -e "\n"

ad_high=$(awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/hightemp.motu | wc -l)
ad_normal=$(awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/normaltemp.motu | wc -l)

if [ $ad_normal > $ad_high ]; then
    echo "Alpha diversity in normal temp ($ad_normal) is greater than in high temp ($ad_high)" 
else
    echo "Alpha diversity in high temp ($ad_high) is greater than in normal temp ($ad_normal)" 
fi
echo -e "\n"


#1.5 Can you detect any algae organisms in the normal-temperature sample? Report if so. If not, why not?

echo 1.5 Can you detect any algae organisms in the normal-temperature sample? Report if so. If not, why not?
echo

#echo motus profile -f metagenomics-hotspring-normaltemp.1.fq.gz -r metagenomics-hotspring-normaltemp.2.fq.gz -o ../results/metagenomics/normaltemp_kingdom.motu -n normaltemp_kingdom -k kindgdom

f_file=$srcpath/metagenomics-hotspring-normaltemp.1.fq.gz
r_file=${f_file/1/2}
# remove prefix
sample=${f_file##$srcpath/metagenomics-hotspring-}
# remove suffix
sample=${sample%%.1.fq.gz}
motus profile -f $f_file -r $r_file -o $dstpath/${sample}_kingdom.motu -n ${sample}_kingdom -k kingdom

echo -e "\n"
echo detected kingdoms in normal temp sample
echo
awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/${sample}_kingdom.motu
echo -e "\n"


motus profile -f $f_file -r $r_file -o $dstpath/${sample}_phylum.motu -n ${sample}_phylum -k phylum
awk 'NR>3 {if ($NF>0 && $1!="-1") {print}}' $dstpath/${sample}_phylum.motu | grep Cyanobacteria



