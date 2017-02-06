#!/bin/bash

G1a=("1000313262-HCV_S1_L001_R1_001.fastq.gz
1000313806-HCV_S2_L001_R1_001.fastq.gz
1000314557-HCV_S3_L001_R1_001.fastq.gz
1000315165-HCV_S4_L001_R1_001.fastq.gz
1000317127-HCV_S5_L001_R1_001.fastq.gz
1000317351-HCV_S6_L001_R1_001.fastq.gz
1000317355-HCV_S7_L001_R1_001.fastq.gz
1000317492-HCV_S8_L001_R1_001.fastq.gz
1000317535-HCV_S9_L001_R1_001.fastq.gz
1000319329-HCV_S10_L001_R1_001.fastq.gz
1000320219-HCV_S12_L001_R1_001.fastq.gz
1000320509-HCV_S11_L001_R1_001.fastq.gz")

G1b=("")

G3a=("1000335439-SSIV-HCV_S14_L001_R1_001.fastq.gz
1000335439-Takara-HCV_S13_L001_R1_001.fastq.gz")

G4d=("")

HIV2=("1000326309-HIV2_S15_L001_R1_001.fastq.gz
pos-control-HIV2_S16_L001_R1_001.fastq.gz")

n=50000 # number of reads (default 50000)
c=10     # minimal coverage (default 3)
i=4     # iterations (default 3)


script_dir=$( dirname "$(readlink -f "$0")" )

for s in $G1a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_1.fasta -n $n -c $c -i $i $s
done

for s in $G1b; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_1b_M1LE.fasta -n $n -c $c -i $i $s
done

for s in $G3a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_3a_CB.fasta -n $n -c $c -i $i $s
done

for s in $G4d; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_4d.fasta -n $n -c $c -i $i $s
done

for s in $HIV2; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-2A_ROD.fasta -n $n -c $c -i $i $s
done

#Rscript ${script_dir}/cov_plot.R
#Rscript ${script_dir}/spreadvcf