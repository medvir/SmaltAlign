#!/bin/bash

HCV-a=("1000331007_S2_L001_R1_001.fastq.gz
1000333464_S1_L001_R1_001.fastq.gz
1000335027_S3_L001_R1_001.fastq.gz
1000335296_S4_L001_R1_001.fastq.gz
1000335797_S5_L001_R1_001.fastq.gz
1000336164_S6_L001_R1_001.fastq.gz
1000336521_S7_L001_R1_001.fastq.gz
1000336905_S12_L001_R1_001.fastq.gz
1000337353_S8_L001_R1_001.fastq.gz
1000338554_S9_L001_R1_001.fastq.gz
1000338555_S10_L001_R1_001.fastq.gz
1000338688_S11_L001_R1_001.fastq.gz
1000340531_S16_L001_R1_001.fastq.gz
1000340830_S13_L001_R1_001.fastq.gz
1000342005_S14_L001_R1_001.fastq.gz
1000342493_S15_L001_R1_001.fastq.gz")

HCV-1b=("")

HCV-1l=("")

HCV-3a=("")

HCV-4d=("")

HIV-1=("1000342682_S18_L001_R1_001.fastq.gz
1000342899_S19_L001_R1_001.fastq.gz
1000343094_S20_L001_R1_001.fastq.gz
1000343104_S21_L001_R1_001.fastq.gz
1000343237_S22_L001_R1_001.fastq.gz
1000343319_S23_L001_R1_001.fastq.gz
1000343328_S24_L001_R1_001.fastq.gz
1000343615_S25_L001_R1_001.fastq.gz
1000343652_S17_L001_R1_001.fastq.gz")

HIV-2=("")


n=50000 # number of reads (default 50000)
i=4     # iterations (default 4)
script_dir=$( dirname "$(readlink -f "$0")" )


for s in $HCV-1a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1a.fasta -n $n -i $i $s
done


for s in $HCV-1b; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1b.fasta -n $n -i $i $s
done


for s in $HCV-1l; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1l.fasta -n $n -i $i $s
done


for s in $HCV-3a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-3a.fasta -n $n -i $i $s
done


for s in $HCV-4d; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-4d.fasta -n $n -i $i $s
done


for s in $HIV-1; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-1.fasta -n $n -i $i $s
done


for s in $HIV-2; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-2.fasta -n $n -i $i $s
done


#Rscript ${script_dir}/cov_plot.R ./
#Rscript ${script_dir}/wts.R ./
