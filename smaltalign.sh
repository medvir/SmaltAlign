#!/bin/bash

# arguments:
# 1. fastq file or directory containing fastq files
# 2. reference in fasta format
# [3. number of reads]


#################################

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
backbone=${script_dir}/pcDNA3_bb.fasta
reference=${script_dir}/HXB2.fasta
reads_limit=100000
expected_length=3500

### arguments
if [ $# == 0 ]; then
	echo
	echo 'enveq.sh [options] ...'
	echo
	echo '-r, --reference        reference (default HXB2.fasta)'
	echo '-b, --backbone         plasmid backbone to be subtracted (default pcDNA3_bb.fasta'
	echo '-n, --reads_limit      limit number of reads (default 100000)'
	echo '-l, --expected_length  expected insert length (default 3500)'
	echo
	exit
fi

while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-r|--reference)
		reference="$2"
		shift # past argument
		;;
		-b|--backbone)
		backbone="$2"
		shift # past argument
		;;
		-n|--reads_limit)
		reads_limit="$2"
		shift # past argument
		;;
		-l|--expected_length)
		expected_length="$2"
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
backbone=$( readlink -f $backbone )
reference=$( readlink -f $reference )

echo
echo -e 'sample_dir \t' $sample_dir
echo -e 'script_dir \t' $script_dir
echo -e 'backbone \t' $backbone
echo -e 'reference \t' $reference
echo -e 'reads_limit \t' $reads_limit
echo -e 'expected_length ' $expected_length
echo







############################




ref=$2

### make list of files to analyse
if [ -d $1 ]; then list=$(ls $1 | grep .fastq); else list=$1; fi

### set n = 1 if third argument not given
if [ -z $3 ]; then n=1; else n=$3; fi

### loop over all files	
for i in $list; do
	
	### unzip fastq file if necessary
	if [[ $i =~ \.gz$ ]]
		then
			gunzip $i
			i=$(echo $i | sed 's/.gz//')
			echo $i
		fi
	
	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/.fastq//')
	
	echo 
	echo $name
	echo ----------------------------------------


	### sample subset of reads
	seqtk sample $i $n > ${name}_reads.fastq	

	### align with smalt to refernce
	smalt index -k 7 -s 2 smalt_index $ref
	smalt map -n 28 -x -f samsoft -o ${name}.sam smalt_index ${name}_reads.fastq
	samtools view -Su ${name}.sam | samtools sort - ${name}
	samtools index ${name}.bam	

	### create consensus
	freebayes -f $ref -p 1 ${name}.bam > ${name}.vcf	
	vcf2fasta -f $ref -p ${name}_ -P 1 ${name}.vcf
	mv ${name}_unknown* ${name}_cons.fasta

	### create vcf with lofreq
	rm ${name}_lofreq.vcf
	lofreq call -f $ref -o ${name}_lofreq.vcf ${name}.bam

	### calculate depth
	samtools depth ${name}.bam | cut -f 2,3 > ${name}.depth

	### remove temporary files
	rm ${name}.sam
	#rm ${name}.vcf
	rm ${name}_reads.fastq

	### run gap_cons.py
	/home/huber.michael/GapCons/gap_cons.py ${name}_cons.fasta ${name}.depth 1

done
