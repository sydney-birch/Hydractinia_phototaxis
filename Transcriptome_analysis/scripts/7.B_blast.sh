#! /bin/bash

#This script will blast your query file to your human fasta and will output a dir called blastout with your result. This particular blast script will only find 
#1 match for each gene in your query file at an evalue of 0.00001. 

#fasta_query=$1

for fasta in ./blastdb/*.fa
do
    echo "BLASTING  $fasta ..."
    blastp -query $1 -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001
done

mkdir blastout

mv ./blastdb/*blastout ./blastout/
