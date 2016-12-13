#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
n=100000
min_cov=3

### arguments
if [ $# == 0 ]; then
	echo
	echo 'Usage: smaltalign.sh -r <reference_file> [options] <fastq_file/directory> '
	echo 'Options:'
	echo '  -r       reference'
	echo '  -n INT   number of reads (default 100,000)'
	echo '  -c INT   minimal coverage (default 3)'
	echo
	exit
fi

while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-r)
		ref="$2"
		shift # past argument
		;;
		-c)
		min_cov="$2"
		shift # past argument
		;;
		-n)
		n="$2"
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
ref=$( readlink -f $ref )
sample_dir=$( readlink -f $sample_dir )

### print arguments
echo -e 'sample_dir: ' $sample_dir
echo -e 'script_dir: ' $script_dir
echo -e 'ref: ' $ref
echo -e 'n: ' $n
echo -e 'min_cov: ' $min_cov

### make list of files to analyse
if [ -d $sample_dir ]; then list=$(ls $sample_dir | grep .fastq); else list=$sample_dir; fi

### loop over list	
for i in $list; do
	
	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/.fastq.gz//'| sed 's/.fastq//')
	echo $name
	echo ----------------------------------------
	
	### unzip fastq file if necessary, sample reads with seqtk
	if [[ $i =~ \.gz$ ]]
		then
			gunzip -c $i > ${name}_unzipped.fastq
			seqtk sample ${name}_unzipped.fastq $n > ${name}_reads.fastq
		else
			seqtk sample $i $n > ${name}_reads.fastq
		fi
		
	### align with smalt to reference
	smalt index -k 7 -s 2 smalt_index $ref
	smalt map -n 28 -x -f samsoft -o ${name}.sam smalt_index ${name}_reads.fastq
	samtools view -Su ${name}.sam | samtools sort - ${name}
	samtools index ${name}.bam	

	### create consensus wiht freebayes
	freebayes -f $ref -p 1 ${name}.bam > ${name}.vcf	
	vcf2fasta -f $ref -p ${name}_ -P 1 ${name}.vcf
	mv ${name}_unknown* ${name}_cons.fasta

	### create vcf with lofreq
	rm ${name}_lofreq.vcf
	lofreq call -f $ref -o ${name}_lofreq.vcf ${name}.bam
	
	### calculate depth, run gap_cons.py
	samtools depth ${name}.bam > ${name}.depth
	$script_dir/gapcons.py ${name}_cons.fasta ${name}.depth $min_cov
	
	### remove temporary files
	rm ${name}.sam
	rm ${name}.vcf
	rm ${name}_reads.fastq
	rm ${name}_unzipped.fastq
	
done