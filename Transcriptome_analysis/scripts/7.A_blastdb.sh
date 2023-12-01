#!/bin/bash

#make a blast database. You will need to make a dir called fastas and inside that dir add your human fasta. 

for i in ./fastas/*.fa
do
    echo "formating for BLAST  $i ..."
    makeblastdb -in $i -parse_seqids -dbtype prot
done

# move blastdbs to new dir
mkdir blastdb
cp ./fastas/* blastdb

# clean up fastas dir
rm ./fastas/*.fa.*
