#!/bin/bash

HCV-1a=("1000331007_S2_L001_R1_001.fastq.gz")

HCV-1b=("")

HCV-1l=("")

HCV-3a=("")

HCV-4d=("")

HIV-1=("1000342682_S18_L001_R1_001.fastq.gz
1000343615_S25_L001_R1_001.fastq.gz
1000343652_S17_L001_R1_001.fastq.gz")

HIV-2=("")


n=200000 # number of reads (default 50000)
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
