#!/bin/bash

# set script directory
#script_dir=$(dirname "$(readlink -f "$0")")
script_dir="/home/schmutz.stefan/Repositories/SmaltAlign/"

# set sample directory
sample_dir=$(dirname "$(readlink -f "$0")")

# set SmaltAlign parameters
n=200000 # number of reads (default 200000)
i=4 # iterations (default 4)

# define Virus and subtype of the sequencing samples
HCV_1a=("")

HCV_1b=("")

HCV_1l=("")

HCV_2c=("")

HCV_3a=("")

HCV_4d=("")

HIV_1=("")

HIV_2=("")

HSV_1=("")

# run SmaltAlign for all defined samples
for s in $HCV_1a; do
    ${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1a.fasta -n $n -i $i $s
done

for s in $HCV_1b; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1b.fasta -n $n -i $i $s
done

for s in $HCV_1l; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-1l.fasta -n $n -i $i $s
done

for s in $HCV_2c; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-2c.fasta -n $n -i $i $s
done

for s in $HCV_3a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-3a.fasta -n $n -i $i $s
done

for s in $HCV_4d; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV-4d.fasta -n $n -i $i $s
done

for s in $HIV_1; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-1.fasta -n $n -i $i $s
done

for s in $HIV_2; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-2.fasta -n $n -i $i $s
done

for s in $HSV_1; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HSV-1.fasta -n $n -i $i $s
done

# run R scripts
Rscript ${script_dir}/cov_plot.R ${sample_dir}
Rscript ${script_dir}/wts.R ${sample_dir}
