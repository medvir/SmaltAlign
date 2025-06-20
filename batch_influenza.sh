#!/bin/bash

set -o pipefail

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: 
batch_influenza.sh [options]
Options:
    [-h or --help]
    [-r or --refdir]
    [-s or --sampledir]
    [-n or --numreads]
    [-i or --iterations]
    [-t or --varthres]
    [-c or --mincov]
    [-o or --outdir]
"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
Options:
    -h, --help:
        Show this help message and exit.
    -r, --refdir:
        Path to the reference directory needed for the analysis.
        The file must contain the references for all segments.
        If it doesn't exist, the best reference for each sample
        and segment will be chosen automatically.
        Default: ./Reference_sequences.
    -s, --sampledir:
        Path to the sample directory.
        It must contain one (or more) FASTQ files.
        Names of the FASTQ files must end with ".fastq.gz".
        Default: current directory.
    -n, --numreads:
        Number of reads used for subsampling.
        Default: 10000000.
    -i, --iterations:
        Number of iterations.
        Default: 4.
    -t, --varthres:
        Minority variant threshold, percentage.
        Default: 15.
    -c, --mincov:
        Minimal coverage required (reads).
        Default: 3.
    -o, --outdir:
        Path to the output directory (it must exist).
        Default: current directory.
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--refdir") set -- "$@" "-r" ;;
                "--sampledir") set -- "$@" "-s" ;;
                "--numreads") set -- "$@" "-n" ;;
                "--iterations") set -- "$@" "-i" ;;
                "--varthres") set -- "$@" "-t" ;;
                "--mincov") set -- "$@" "-c" ;;
                "--outdir") set -- "$@" "-o" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
outdir="./"; n=10000000; i=4
ref_dir="./Reference_sequences"
sample_dir="./"
varthres=15; mincov=3

# Define all parameters
while getopts 'r::s::n::i::t::c::o::dh' flag; do
        case "${flag}" in
                r) ref_dir=${OPTARG} ;;
                s) sample_dir=${OPTARG} ;;
                n) n=${OPTARG} ;;
                i) i=${OPTARG} ;;
                t) varthres=${OPTARG} ;;
                c) mincov=${OPTARG} ;;
				o) outdir=${OPTARG} ;;
                h) print_help
                   exit 1;;
                *) print_usage
                    exit 1;;
        esac
done

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )

### convert relative to absolute path
sample_dir=$( readlink -f $sample_dir )
outdir=$( readlink -f $outdir )

for filename in $sample_dir/*.fastq.gz; do

    name=$(echo $(basename "$filename") | cut -f 1 -d '.')

    for segment in {1..8}; do
        mkdir -p ${outdir}/${name}/segment${segment}
        mkdir ${outdir}/${name}/segment${segment}/cons
    done

    # choose best reference for each sample and segment
    if [ -d "$ref_dir" ]; then
        ref_dir=$( readlink -f $ref_dir )
    else
        mkdir ${outdir}/${name}/Reference_sequences
        python $script_dir/select_ref_segments.py $filename ${outdir}/${name}/Reference_sequences ${script_dir}
        ref_dir=$( readlink -f ${outdir}/${name}/Reference_sequences )
    fi

    for segment in {1..8}; do
    (   # run smaltalign.sh with the previously selected best reference
        ${script_dir}/smaltalign.sh \
        -r $ref_dir/*segment*$segment*fasta \
        -n $n \
        -i $i \
        -t $varthres \
        -c $mincov \
        -o ${outdir}/${name}/segment${segment} $filename
        
        cp ${outdir}/${name}/segment${segment}/*WTS.fasta ${outdir}/${name}/segment${segment}/cons/ )
    done
    wait
done
