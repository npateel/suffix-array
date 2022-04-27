#!/bin/bash

SIZES="10 100 1000 10000 100000 1000000 10000000"

#wget -nc  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

#pigz -d GRCh38_latest_genomic.fna.gz
#sed 's/N//g' GRCh38_latest_genomic.fna | sed '/^$/d' | tr [:lower:] [:upper:] > human.fna 
#sed 's/N//g' ecoli.fa | sed '/^$/d' | tr [:lower:] [:upper:] > human.fna 


for number in `echo $SIZES `
do
	head -n $number ecoli.fa > "samples/ecoli-$number.fa"
done
