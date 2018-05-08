#!/bin/bash

# set directory where smaltalign.sh is located
script_dir="/analyses/Diagnostics/Repositories/SmaltAlign/"

# set sample directory to the directory where this file is
sample_dir=$(dirname "$(readlink -f "$0")")

# set SmaltAlign parameters
n=200000 # number of reads (default 200000)
i=4 # iterations (default 4)

source activate SmaltAlign

for filename in $sample_dir/*.fastq.gz; do

    name=$(echo "$filename" | cut -f 1 -d '.')

    for segment in {1..8}; do
        mkdir -p $name/segment$segment
    done

    cd $name

    # choose best reference for each sample and segment
    $script_dir/select_ref.py $filename

    for segment in {1..8}; do
    (   cd segment$segment
        # run smaltalign.sh with the previously selected best reference
        ${script_dir}/smaltalign.sh \
        -r ../segment-$segment.fasta \
        -n $n \
        -i $i $filename

        Rscript ${script_dir}/cov_plot.R ./
        source deactivate
        Rscript ${script_dir}/wts.R ./
        source activate SmaltAlign
        cd ../ ) &
    done
    wait

    cd ../
done
