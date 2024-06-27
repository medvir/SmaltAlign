#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
n_reads=200000
iterations=4

### arguments
if [ $# == 0 ]; then
	echo
	echo 'Usage: smaltalign.sh -r <reference_file> [options] <fastq_file/directory> '
	echo 'Options:'
	echo '  -r       reference'
	echo '  -n INT   number of reads (default 200000)'
	echo '  -i INT   iterations (default 4)'
	echo
	exit
fi

while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-r)
		all_ref="$2"
		shift # past argument
		;;
		-n)
		n_reads="$2"
		shift # past argument
		;;
		-i)
		iterations="$2"
		shift # past argument
		;;
		*)
		# unknown option
		;;
	esac
	shift # past argument or value
done

if [[ -n $1 ]]; then
    sample_dir=$1
fi

### convert relative to absolute path
all_ref=$( readlink -f $all_ref )
sample_dir=$( readlink -f $sample_dir )

### print arguments
echo -e 'sample_dir: ' $sample_dir
echo -e 'script_dir: ' $script_dir
echo -e 'ref: ' $all_ref
echo -e 'n_reads: ' $n_reads
echo -e 'iterations: ' $iterations




### loop over list of files to analyse
if [ -d $sample_dir ]; then list=$(ls $sample_dir | grep .fastq); else list=$sample_dir; fi
for i in $list; do

	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/.fastq.gz//'| sed 's/.fastq//')
    
    
	
	### sample reads with seqtk
	seqtk sample $i $n_reads > ${name}_reads.fastq 
	n_sample=$(wc -l ${name}_reads.fastq | cut -f 1 -d " ")
	#n_sample=$(($n_sample / 4))
	
	# select the most probable reference
	echo $all_ref
	echo ${name}_reads.fastq
	python $script_dir/select_reference.py -f ${name}_reads.fastq -r $all_ref -s 1000
	mv ./reference_freq.csv ./${name}_references_freq.csv
	mv ./chosen_reference.fasta ./${name}_chosen_reference.fasta
	ref=${name}_chosen_reference.fasta
	echo $ref

	### de novo aligment
	echo
	echo sample $name, $n_sample reads, de-novo alignment
	echo "*******************************************************************************"
	velveth ${name} 29 -fastq ${name}_reads.fastq
	velvetg ${name} -min_contig_lgth 200

	### fake fastq of contigs and cat to reads in triplicate
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${name}/contigs.fa | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+"$1"\n"gensub(/./, "I", "g", $2)}' > ${name}/contigs.fastq
	cat ${name}_reads.fastq ${name}/contigs.fastq ${name}/contigs.fastq ${name}/contigs.fastq > ${name}_reads_contigs.fasta

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
		smalt index -k 7 -s 2 ${name}_${it}_smalt_index $ref
		samtools faidx $ref

		if [ "$it" -eq 1 ]; then
			smalt map -n 24 -x -y 0.5 -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads_contigs.fasta
		else
			smalt map -n 24 -x -y 0.5 -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads.fastq
		fi

		samtools view -Su ${name}_${it}.sam | samtools sort -o ${name}_${it}.bam
		samtools index ${name}_${it}.bam

		### create consensus with freebayes
		freebayes -f $ref -p 1 ${name}_${it}.bam > ${name}_${it}.vcf
		muts=$(grep -c -v "^#" ${name}_${it}.vcf)
        if [ "$muts" -eq "0" ]; then
            echo WARNING: vcf file empty
            cat $ref > ${name}_${it}_cons.fasta
        else
            vcf2fasta -f $ref -p ${name}_${it}_ -P 1 ${name}_${it}.vcf
            mv ${name}_${it}_unknown* ${name}_${it}_cons.fasta
        fi

		### create vcf with lofreq
		rm -f ${name}_${it}_lofreq.vcf

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
        lofreq call-parallel --pp-threads $cores -f $ref -o ${name}_${it}_lofreq.vcf ${name}_${it}.bam

		### calculate depth
		samtools depth -d 1000000 ${name}_${it}.bam > ${name}_${it}.depth

		### remove temporary files of this iteration
		rm -f ${name}_${it}.sam
		rm -f ${name}_${it}.vcf
		rm -f ${name}_${it}_smalt_index.*
		rm -f ${name}_*_cons.fasta.fai

		### new ref for next iteration
		ref=${name}_${it}_cons.fasta
		((it+=1))
    done

	### remove temporary files
	rm -f ${name}_reads.fastq
	rm -f ${name}_reads_contigs.fasta
	rm -rf ${name}
done
