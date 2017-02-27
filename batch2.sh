#!/bin/bash
	
n=10000 # number of reads (default 50000)
i=4     # iterations (default 4)

G1a=("BE0664-02_S30_L001_R1_001.fastq.gz
BE0665-03_S31_L001_R1_001.fastq.gz
BE2378-02_S33_L001_R1_001.fastq.gz
BE2899-01_S26_L001_R1_001.fastq.gz
BE7586-02_S27_L001_R1_001.fastq.gz
BE7587-01_S28_L001_R1_001.fastq.gz
BE8666-01_S29_L001_R1_001.fastq.gz
BE9308-02_S32_L001_R1_001.fastq.gz
KoNeg_S35_L001_R1_001.fastq.gz
LowKoPos_S34_L001_R1_001.fastq.gz")

script_dir=$( dirname "$(readlink -f "$0")" )

for s in $G1a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/Cucurbita_pepo.fasta -n $n -i $i $s
done