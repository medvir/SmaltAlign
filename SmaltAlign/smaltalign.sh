#!/bin/bash

set -o pipefail

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: 
smaltalign.sh -r <reference_file> [options] <fastq_file/directory>
Options:
	[-h or --help]
	[-n or --numreads]
	[-i or --iterations]
	[-t or --varthres]
	[-c or --mincov]
	[-o or --outdir]
	[-d or --indels]
	[-m or --mergecov]
	[-a or --alltogether]
	[-p or --paired]
"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
Options:
	-h, --help:
		Show this help message and exit.
	-n, --numreads:
		Number of reads used for subsampling.
		Default: 200000.
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
	-d, --indels:
		If required, extract indels as well.
		Default: don't perform this step.
	-m, --mergecov:
		If required, generate an additional covplot merging all results.
		Default: don't perform this step.
	-a, --alltogether:
		Merge all input FASTQ files into a single combined run.
		Outputs include a single consensus sequence and a single coverage plot.
		Default: process files individually.
	-p, --paired:
		Enable paired-end reads processing. It automatically detects 
		separate files (e.g., _R1 and _R2) or an interleaved file.
		Default: single-end mode.
Required arguments:
	-r, --reference:
		Path to the reference file needed for the analysis.
		The file can contain one reference or several.
		If it has several references, the best one is automatically chosen.
	fastq_file/directory:
		Path to the FASTQ file to analyse.
		Note that no flag is required for it.
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--reference") set -- "$@" "-r" ;;
                "--numreads") set -- "$@" "-n" ;;
                "--iterations") set -- "$@" "-i" ;;
				"--varthres") set -- "$@" "-t" ;;
                "--mincov") set -- "$@" "-c" ;;
                "--outdir") set -- "$@" "-o" ;;
				"--indels") set -- "$@" "-d" ;;
				"--mergecov") set -- "$@" "-m" ;;
				"--alltogether") set -- "$@" "-a" ;;
				"--paired") set -- "$@" "-p" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
outdir_final="./"; n_reads=200000; iterations=4; indels=0
varthres=15; mincov=3; mergecov=0; paired=0; alltogether=0

# Define all parameters
while getopts 'r:n::i::t::c::o::dmaph' flag; do
        case "${flag}" in
                r) reference=${OPTARG} ;;
                n) n_reads=${OPTARG} ;;
                i) iterations=${OPTARG} ;;
				t) varthres=${OPTARG} ;;
                c) mincov=${OPTARG} ;;
				o) outdir_final=${OPTARG} ;;
                d) indels=${OPTARG} ;;
				m) mergecov=${OPTARG} ;;
				a) alltogether=1 ;;
				p) paired=1 ;;
                h) print_help
                   exit 1;;
                *) print_usage
                    exit 1;;
        esac
done
shift $(expr $OPTIND - 1 )
if [[ -n $1 ]]; then
    sample_dir=$1
else
	echo "Missing FASTQ file!"
	exit 1
fi

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )

### convert relative to absolute path
reference=$( readlink -f $reference )
sample_dir=$( readlink -f $sample_dir )
outdir_final=$( readlink -f $outdir_final )

### print arguments
echo -e 'sample_dir: ' $sample_dir
echo -e 'script_dir: ' $script_dir
echo -e 'reference: ' $reference
echo -e 'n_reads: ' $n_reads
echo -e 'iterations: ' $iterations
echo -e 'varthres: ' $varthres
echo -e 'mincov: ' $mincov
echo -e 'outdir: ' $outdir_final
echo -e 'indels: ' $indels
echo -e 'mergecov: ' $mergecov
echo -e 'alltogether: ' $alltogether
echo -e 'paired_end: ' $paired
# Specify temporal location for all intermediate files
rand_id=$(tr -dc 'a-zA-Z' < <(date +%s%N | sha256sum) | head -c5)
outdir="/tmp/SmaltAlignOutputs_${rand_id}/"
rm -rf ${outdir}
mkdir ${outdir}

### loop over list of files to analyse
if [[ $alltogether == 1 ]]; then
	if [ -d $sample_dir ]; then
        echo " ALL-TOGETHER MODE ENABLED: Merging all input files..."
        merged_dir="${outdir}/MergedInput"
        mkdir -p "$merged_dir"

		# Initialize a variable to count how many samples we merge
		merge_count=0

        if [[ $paired == 1 ]]; then
			# 1. Check if separated _R1 and _R2files exist
            r1_files=$(ls $sample_dir | grep "\.fastq" | grep "_R1" || true)
			r2_files=$(ls $sample_dir | grep "\.fastq" | grep "_R2" || true)

			if [[ -n "$r1_files" && -n "$r2_files" ]]; then
                echo "Detected separated R1 and R2 files."

				merge_count=$(echo $r1_files | wc -w)
            
				echo "Merging all R1 files..."
				for f in $(ls $sample_dir | grep "\.fastq" | grep "_R1"); do
					if [[ $f == *.gz ]]; then gunzip -c "${sample_dir}/$f"; else cat "${sample_dir}/$f"; fi
				done > "${merged_dir}/MergedAll_R1.fastq"

				echo "Merging all R2 files..."
				for f in $(ls $sample_dir | grep "\.fastq" | grep "_R2"); do
					if [[ $f == *.gz ]]; then gunzip -c "${sample_dir}/$f"; else cat "${sample_dir}/$f"; fi
				done > "${merged_dir}/MergedAll_R2.fastq"
			else
				echo "No separated R1 and R2 files detected. Assuming interleaved files."
				echo "Merging all interleaved files..."
				all_fastqs=$(ls $sample_dir | grep "\.fastq" || true)
				merge_count=$(echo $all_fastqs | wc -w)
				for f in $all_fastqs; do
					if [[ $f == *.gz ]]; then gunzip -c "${sample_dir}/$f"; else cat "${sample_dir}/$f"; fi
				done > "${merged_dir}/MergedAll.fastq"
			fi
        else
            echo "Merging all single-end files..."
			all_fastqs=$(ls $sample_dir | grep "\.fastq" || true)
			merge_count=$(echo $all_fastqs | wc -w)
            for f in $all_fastqs; do
                if [[ $f == *.gz ]]; then gunzip -c "${sample_dir}/$f"; else cat "${sample_dir}/$f"; fi
            done > "${merged_dir}/MergedAll.fastq"
        fi

		# Adjust the number of randomly sampled reads to maintain the same depth across the merged samples
		if [[ $n_reads != "all" ]]; then
            old_reads=$n_reads
            n_reads=$(( n_reads * merge_count ))
            echo "Scaled downsampling (-n) from $old_reads to $n_reads to maintain depth across $merge_count combined samples."
        fi
        
        # Override the sample directory to consider only the merged files
        sample_dir=$merged_dir
    else
        echo "WARNING: -a/--alltogether was selected, but a single file was provided instead of a directory. Ignoring merge."
    fi
fi

if [ -d $sample_dir ]; then 
	if [[ $paired == 1 ]]; then
		list=$(ls $sample_dir | grep "\.fastq" | grep -v "_R2" | sed -e "s#^#${sample_dir}/#")
	else
		list=$(ls $sample_dir | grep .fastq | sed -e "s#^#${sample_dir}/#")
	fi 
else 
	list=${sample_dir}
fi
num_files=$(echo $list | wc -w)

# Calculate iterations that will be later removed (intermediate files)
minus_one=$((iterations - 1))

for i in $list; do
	
	ref=$reference
	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/_R[12].*//' | sed 's/.fastq.gz//'| sed 's/.fastq//')
	if [[ $mergecov != 0 ]]; then
		all_names="${all_names}${name} "
	fi

	if [ $num_files -gt 1 ]; then
		mkdir ${outdir}/${name}
		new_outdir="${outdir}/${name}/"
	else
		new_outdir="${outdir}"
	fi

	### File Setup and Subsampling
    if [[ $paired == 1 ]]; then
        # 1. Determine if separate or interleaved
        r2_file=""
        if [[ $i == *"_R1"* ]]; then
            r2_file="${i/_R1/_R2}"
        fi

        if [[ -n "$r2_file" && -f "$r2_file" ]]; then
            echo "Paired-end reads detected: Separate files"
            r1_input=$i
            r2_input=$r2_file
            is_interleaved=0
        else
            echo "Paired-end reads detected: Interleaved file"
            # De-interleave temporary files using awk to properly sample pairs
            if [[ $i == *.gz ]]; then
                gunzip -c $i | awk 'NR%8==1||NR%8==2||NR%8==3||NR%8==4' > ${new_outdir}/${name}_R1_full.fastq
                gunzip -c $i | awk 'NR%8==5||NR%8==6||NR%8==7||NR%8==0' > ${new_outdir}/${name}_R2_full.fastq
            else
                awk 'NR%8==1||NR%8==2||NR%8==3||NR%8==4' $i > ${new_outdir}/${name}_R1_full.fastq
                awk 'NR%8==5||NR%8==6||NR%8==7||NR%8==0' $i > ${new_outdir}/${name}_R2_full.fastq
            fi
            r1_input="${new_outdir}/${name}_R1_full.fastq"
            r2_input="${new_outdir}/${name}_R2_full.fastq"
            is_interleaved=1
        fi

		### sample reads with seqtk (subsample with identical random seeds to preserve pairing)
		seed=$RANDOM
		if [[ $n_reads == "all" ]]; then 
			if [[ $r1_input == *.gz ]]; then 
				gunzip -c $r1_input > ${new_outdir}/${name}_R1_reads.fastq
			else 
				cp $r1_input ${new_outdir}/${name}_R1_reads.fastq
			fi
			if [[ $r2_input == *.gz ]]; then 
				gunzip -c $r2_input > ${new_outdir}/${name}_R2_reads.fastq
			else 
				cp $r2_input ${new_outdir}/${name}_R2_reads.fastq
			fi
			n_sample="all"
		else
			seqtk sample -s $seed $r1_input $n_reads > ${new_outdir}/${name}_R1_reads.fastq
            seqtk sample -s $seed $r2_input $n_reads > ${new_outdir}/${name}_R2_reads.fastq
            n_sample=$(wc -l ${new_outdir}/${name}_R1_reads.fastq | cut -f 1 -d " ")
            n_sample=$(($n_sample / 4))
		fi

		if [[ $is_interleaved == 1 ]]; then
            rm -f ${new_outdir}/${name}_R1_full.fastq ${new_outdir}/${name}_R2_full.fastq
        fi
	else
        # Single-end handling
        if [[ $n_reads == "all" ]]; then 
            if [[ $i == *.gz ]]; then 
				gunzip -c $i > ${new_outdir}/${name}_reads.fastq
			else 
				cp $i ${new_outdir}/${name}_reads.fastq
			fi
            n_sample="all"
        else
            seqtk sample $i $n_reads > ${new_outdir}/${name}_reads.fastq
            n_sample=$(wc -l ${new_outdir}/${name}_reads.fastq | cut -f 1 -d " ")
            n_sample=$(($n_sample / 4))
        fi
    fi

	# select the most probable reference
	n_refs=$(grep "^>" $ref | wc -l)
	if [ "$n_refs" -gt 1 ]; then
		if [[ $paired == 1 ]]; then
            python $script_dir/select_ref_whole.py -f ${new_outdir}/${name}_R1_reads.fastq -r $ref -s 1000 -o $new_outdir
        else
			python $script_dir/select_ref_whole.py -f ${new_outdir}/${name}_reads.fastq -r $ref -s 1000 -o $new_outdir
		fi
		mv ${new_outdir}/reference_freq.csv ${new_outdir}/${name}_references_freq.csv
		mv ${new_outdir}/chosen_reference.fasta ${new_outdir}/${name}_chosen_reference.fasta
		ref=${new_outdir}/${name}_chosen_reference.fasta
	fi
	seqkit replace -p "/|:" -r '_' $ref > ${new_outdir}/${name}_chosen_mod_reference.fasta
	mv ${new_outdir}/${name}_chosen_mod_reference.fasta ${new_outdir}/${name}_chosen_reference.fasta
	ref=${new_outdir}/${name}_chosen_reference.fasta
	echo "Using reference: $ref"

	### de novo aligment
	echo
	echo sample $name, $n_sample reads, de-novo alignment
	echo "*******************************************************************************"

	if [[ $paired == 1 ]]; then
        velveth ${new_outdir}/${name} 29 -fastq -shortPaired -separate ${new_outdir}/${name}_R1_reads.fastq ${new_outdir}/${name}_R2_reads.fastq
    else
		velveth ${new_outdir}/${name} 29 -fastq ${new_outdir}/${name}_reads.fastq
	fi
	velvetg ${new_outdir}/${name} -min_contig_lgth 200

	### fake fastq of contigs and cat to reads in triplicate
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${new_outdir}/${name}/contigs.fa | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+"$1"\n"gensub(/./, "I", "g", $2)}' > ${new_outdir}/${name}/contigs.fastq
	if [[ $paired == 1 ]]; then
        cat ${new_outdir}/${name}_R1_reads.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq > ${new_outdir}/${name}_R1_reads_contigs.fasta
        cat ${new_outdir}/${name}_R2_reads.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq > ${new_outdir}/${name}_R2_reads_contigs.fasta
    else
		cat ${new_outdir}/${name}_reads.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq ${new_outdir}/${name}/contigs.fastq > ${new_outdir}/${name}_reads_contigs.fasta
	fi

	it=1
	while [ "$it" -le "$iterations" ]; do
		echo
		echo sample $name, $n_sample reads, iteration $it
		if [ -s $ref ]; then
			echo reference $ref
			echo "*******************************************************************************"
		else
			echo ERROR: No reference
			exit
		fi

		### align with smalt to reference, reads either ${name}_reads.fastq or ${name}_reads_contigs.fasta
		smalt index -k 7 -s 2 ${new_outdir}/${name}_${it}_smalt_index $ref
		samtools faidx $ref

		if [ "$it" -eq 1 ] && [ "$iterations" -gt 1 ]; then
			if [[ $paired == 1 ]]; then
                smalt map -n 24 -x -y 0.5 -f samsoft -o ${new_outdir}/${name}_${it}.sam ${new_outdir}/${name}_${it}_smalt_index ${new_outdir}/${name}_R1_reads_contigs.fasta ${new_outdir}/${name}_R2_reads_contigs.fasta
            else
				smalt map -n 24 -x -y 0.5 -f samsoft -o ${new_outdir}/${name}_${it}.sam ${new_outdir}/${name}_${it}_smalt_index ${new_outdir}/${name}_reads_contigs.fasta
			fi
		else
			if [[ $paired == 1 ]]; then
                smalt map -n 24 -x -y 0.5 -f samsoft -o ${new_outdir}/${name}_${it}.sam ${new_outdir}/${name}_${it}_smalt_index ${new_outdir}/${name}_R1_reads.fastq ${new_outdir}/${name}_R2_reads.fastq
            else
				smalt map -n 24 -x -y 0.5 -f samsoft -o ${new_outdir}/${name}_${it}.sam ${new_outdir}/${name}_${it}_smalt_index ${new_outdir}/${name}_reads.fastq
			fi
		fi

		samtools view -Su ${new_outdir}/${name}_${it}.sam | samtools sort -o ${new_outdir}/${name}_${it}.bam
		samtools index ${new_outdir}/${name}_${it}.bam

		### create consensus with freebayes
		freebayes -f $ref -p 1 ${new_outdir}/${name}_${it}.bam > ${new_outdir}/${name}_${it}.vcf
		muts=$(grep -c -v "^#" ${new_outdir}/${name}_${it}.vcf)
		if [ "$muts" -eq "0" ]; then
            echo WARNING: vcf file empty
            cat $ref > ${new_outdir}/${name}_${it}_cons.fasta
        else
            vcf2fasta -f $ref -p ${new_outdir}/${name}_${it}_ -P 1 ${new_outdir}/${name}_${it}.vcf
			seqkit replace -p ":" -r '_' ${new_outdir}/${name}_${it}_unknown* > ${new_outdir}/${name}_${it}_cons.fasta
            rm ${new_outdir}/${name}_${it}_unknown*
        fi

		### create vcf with lofreq
		rm -f ${new_outdir}/${name}_${it}_lofreq.vcf

        if [ "$(uname)" == "Darwin" ]; then
			cores=$(sysctl -n hw.ncpu) ### Mac OS X platform
        elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
			cores=$(nproc) ### GNU/Linux platform
        else
            cores=2 ### unkonwn/other platforms
        fi
        cores=$(expr $cores / 2) ### only run on half the cores
		cores="${cores/0/1}"
        echo =-=-= lofreq with $cores cores =-=-=
		lofreq call-parallel --pp-threads $cores -f $ref -o ${new_outdir}/${name}_${it}_lofreq.vcf ${new_outdir}/${name}_${it}.bam

		if [[ $indels != 0 ]]; then
			echo "Extracting indels"
			samtools sort ${new_outdir}/${name}_${it}.bam > ${new_outdir}/${name}_${it}_sorted.bam
			lofreq indelqual --dindel -f $ref ${new_outdir}/${name}_${it}_sorted.bam  -o ${new_outdir}/${name}_${it}_indel.bam
			samtools index -b ${new_outdir}/${name}_${it}_indel.bam
			lofreq call-parallel --pp-threads $cores --call-indels -f $ref -o ${new_outdir}/${name}_${it}_lofreq_indel.vcf ${new_outdir}/${name}_${it}_indel.bam

			if [ $(grep -v "^#" ${new_outdir}/${name}_${it}_lofreq_indel.vcf | wc -l) -gt 0 ] ; then
				lofreq filter --only-indels -a 0.15 -v 10 --indelqual-thresh 20 -i ${new_outdir}/${name}_${it}_lofreq_indel.vcf -o ${new_outdir}/${name}_${it}_lofreq_indel_temp_hq.vcf
				echo lofreq filter is done ${new_outdir}/${name}_${it}_lofreq_indel_temp_hq.vcf 
				if [ $(grep -v "^#" ${new_outdir}/${name}_${it}_lofreq_indel_temp_hq.vcf | wc -l) -gt 0 ] ; then
					cp ${new_outdir}/${name}_${it}_lofreq_indel_temp_hq.vcf ${new_outdir}/${name}_lofreq_indel_hq.vcf
					echo ${new_outdir}/lofreq_indel_hq.vcf file is created. 
				fi
			fi
		fi

		### calculate depth
		samtools depth -d 1000000 ${new_outdir}/${name}_${it}.bam > ${new_outdir}/${name}_${it}.depth

		### remove temporary files of this iteration
		rm -f ${new_outdir}/${name}_${it}.sam
		rm -f ${new_outdir}/${name}_${it}.vcf
		rm -f ${new_outdir}/${name}_${it}_smalt_index.*
		rm -f ${new_outdir}/${name}_*_cons.fasta.fai

		### new ref for next iteration
		ref=${new_outdir}/${name}_${it}_cons.fasta
		((it+=1))
	done

	### remove temporary files
	if [[ $paired == 1 ]]; then
        rm -f ${new_outdir}/${name}_R1_reads.fastq ${new_outdir}/${name}_R2_reads.fastq
        rm -f ${new_outdir}/${name}_R1_reads_contigs.fasta ${new_outdir}/${name}_R2_reads_contigs.fasta
    else
		rm -f ${new_outdir}/${name}_reads.fastq
		rm -f ${new_outdir}/${name}_reads_contigs.fasta
	fi
	rm -rf ${new_outdir}/${name}
	
	Rscript ${script_dir}/cov_plot.R ${new_outdir}
    Rscript ${script_dir}/wts.R ${new_outdir} $varthres $mincov

	it=1
	while [ "$it" -le "$minus_one" ]; do
		rm ${new_outdir}/${name}_${it}*
		((it+=1))
	done
done

if [[ $mergecov != 0 ]]; then
	Rscript ${script_dir}/cov_plot.R ${outdir} "${all_names}"
fi

# Move all outputs to the desired output directory
mv ${outdir}/* ${outdir_final}/
rm -rf ${outdir}