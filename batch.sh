#!/bin/bash

G1a=("1000310031_S1_L001_R1_001.fastq.gz
1000310742_S2_L001_R1_001.fastq.gz
1000311130_S3_L001_R1_001.fastq.gz
1000311660_S4_L001_R1_001.fastq.gz
1000312141_S5_L001_R1_001.fastq.gz
1000312484_S6_L001_R1_001.fastq.gz
1000312761_S7_L001_R1_001.fastq.gz")
G1b=("1000312788_S8_L001_R1_001.fastq.gz")
G3a=("1000331002_S9_L001_R1_001.fastq.gz")
G4d=("1000332000_S10_L001_R1_001.fastq.gz")

n=75000
c=10
i=3

for s in $G1a; do
	~/Repositories/SmaltAlign/smaltalign.sh -r ../../references/HCV_1.fasta -n $n -c $c -i $i $s
done

for s in $G1b; do
	~/Repositories/SmaltAlign/smaltalign.sh -r ../../references/HCV_1b_M1LE.fasta -n $n -c $c -i $i $s
done

for s in $G3a; do
	~/Repositories/SmaltAlign/smaltalign.sh -r ../../references/HCV_3a_CB.fasta -n $n -c $c -i $i $s
done

for s in $G4d; do
	~/Repositories/SmaltAlign/smaltalign.sh -r ../../references/HCV_4d.fasta -n $n -c $c -i $i $s
done