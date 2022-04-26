#!/bin/bash

SIZES="10 100 1000 10000 100000 1000000 10000000"

wget -nc -q https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

pigz -d GRCh38_latest_genomic.fna.gz

for number in `echo $SIZES `
do
	head -n $number GRCh38_latest_genomic.fna > "samples/human-$number.fa"
done
