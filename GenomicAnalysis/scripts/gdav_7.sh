#!/bin/bash
# script to get the results from final_project


# 1st Create a log file with all stdout and stderr from the script keeping 
# sourcing bash_profile (needed for conda initialization)
# If sourced from main script initialization is done. No need to source .bash_profile
source ~/scripts/init_gdav $1

conda activate base

# 7. Phylogenetic analysis

allproteomepath=~/final_project
allproteomefile=$allproteomepath/all_reference_proteomes.faa

srcpath=~/results/genome_analysis/functional_analysis
dstpath=~/results/genome_analysis/phylogenetic_analysis
imgspath=$dstpath/trees
mkdir -p $imgspath

# needed for ete3 --image to work
# Generating an image should work according to http://etetoolkit.org/, but raise the error "qt.qpa.xcb: could not connect to display" 
# we don't have an X11 display. Resolved with QT_QPA_PLATFORM="vnc" https://github.com/matplotlib/matplotlib/issues/19360
export QT_QPA_PLATFORM="vnc"
 
echo 7. Phylogenetic analysis

makeblastdb -dbtype prot -in $allproteomefile -out $dstpath/all_pro_blast_db

for file in $srcpath/*.faa
do
    # remove prefix
    outfile=${file##$srcpath/}
    # remove suffix
    outfile=${outfile%%.faa}
    outhomologs=$dstpath/${outfile}_homologs
    blastp -task blastp -query "$file" -db $dstpath/all_pro_blast_db -evalue 0.001 -outfmt 6 -out $outhomologs.blastout

    # Create a FASTA file with all the sequences of selected hits, Tip2: Include the query protein in the same FASTA file
    # First insert query protein
    cp $file $outhomologs.faa
    # Now add all sequences of selected hits
    while read -r genid || [[ -n "${genid}" ]]; do
        sed -n '/'$genid'/,/^>/p' $allproteomefile | sed '${/^>/d}' >> $outhomologs.faa
    done <<< "$(cut -f2  $outhomologs.blastout | uniq)"

    # Build a phylogenetic tree out of the FASTA file
    mafft $outhomologs.faa > $outhomologs.alg
    iqtree -s $outhomologs.alg -m LG --fast
    
   	ete3 view -t $outhomologs.alg.treefile --image "${imgspath}/${outfile}.png" -mc 

    # Generate tree in txt format 
    #ete3 view --text -t $outhomologs.alg.treefile > "${imgspath}/${outfile}_tree.txt"
    #for format in png pdf svg  
    #do
    #	ete3 view -t $outhomologs.alg.treefile --image "${imgspath}/${outfile}_r.$format" -mr 
    #	ete3 view -t $outhomologs.alg.treefile --image "${imgspath}/${outfile}_c.$format" -mc 
    #done 

done


