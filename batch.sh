#!/bin/bash

# set directory where smaltalign.sh is located
script_dir="/analyses/Diagnostics/Repositories/SmaltAlign/"

# set sample directory to the directory where this file is
sample_dir=$(dirname "$(readlink -f "$0")")

# set SmaltAlign parameters
n=200000 # number of reads (default 200000)
i=4 # iterations (default 4)

# define Virus and subtype of the sequencing samples
HCV_1a=("")

HCV_1b=("")

HCV_1h=("")

HCV_1l=("")

HCV_2c=("")

HCV_3a=("")

HCV_3b=("")

HCV_4d=("")

HCV_4f=("")

HIV_1=("")

HIV_2=("")

HSV_1=("")

# define an array containing all Viruses
viruses=(HCV_1a HCV_1b HCV_1h HCV_1l
         HCV_2c
         HCV_3a HCV_3b
         HCV_4d HCV_4f
         HIV_1 HIV_2
         HSV_1)

# run SmaltAlign for all defined samples
for virus in "${viruses[@]}"
do
    for sample in ${!virus}
    do
        if [ "${#sample}" -eq "0" ]
        then
            continue
        fi
        ${script_dir}/smaltalign.sh \
        -r ${script_dir}/References/${virus}.fasta \
        -n $n \
        -i $i $sample
    done
done

# run R scripts with sample_dir as variable input
Rscript ${script_dir}/cov_plot.R ${sample_dir}
Rscript ${script_dir}/wts.R ${sample_dir}
