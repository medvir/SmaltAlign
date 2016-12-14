#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
n_reads=100000
min_cov=3
iterations=3

### arguments
if [ $# == 0 ]; then
	echo
	echo 'Usage: smaltalign.sh -r <reference_file> [options] <fastq_file/directory> '
	echo 'Options:'
	echo '  -r       reference'
	echo '  -n INT   number of reads (default 100,000)'
	echo '  -c INT   minimal coverage (default 3)'
	echo '  -i INT   iterations (default 2)'
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
ref=$( readlink -f $ref )
sample_dir=$( readlink -f $sample_dir )

### print arguments
echo -e 'sample_dir: ' $sample_dir
echo -e 'script_dir: ' $script_dir
echo -e 'ref: ' $ref
echo -e 'n_reads: ' $n_reads
echo -e 'min_cov: ' $min_cov
echo -e 'iterations: ' $iterations

### make list of files to analyse
if [ -d $sample_dir ]; then list=$(ls $sample_dir | grep .fastq); else list=$sample_dir; fi

### loop over list	
for i in $list; do
	
	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/.fastq.gz//'| sed 's/.fastq//')
	
	### unzip fastq file if necessary, sample reads with seqtk
	if [[ $i =~ \.gz$ ]]
		then
			gunzip -c $i | seqtk sample - $n_reads > ${name}_reads.fastq
		else
			seqtk sample $i $n_reads > ${name}_reads.fastq
		fi
		
	### de novo aligment
	#velveth ${name} 29 -fastq ${name}_reads.fastq
	#velvetg ${name}
	#seqtk seq -A ${name}_reads.fastq > ${name}_reads.fasta
	#cat ${name}_reads.fasta ${name}/contigs.fa ${name}/contigs.fa ${name}/contigs.fa > ${name}_reads_contigs.fasta
		
	it=1
	while [ "$it" -le "$iterations" ]; do		
		echo
		echo sample $name iteration $it
		echo "***************************"
		
		### align with smalt to reference
		smalt index -k 7 -s 2 smalt_index $ref
		smalt map -n 28 -x -f samsoft -o ${name}_${it}.sam smalt_index ${name}_reads.fastq #### ${name}_reads.fastq or ${name}_reads_contigs.fasta
		samtools view -Su ${name}_${it}.sam | samtools sort - ${name}_${it}
		samtools index ${name}_${it}.bam	

		### create consensus wiht freebayes
		freebayes -f $ref -p 1 ${name}_${it}.bam > ${name}_${it}.vcf	
		vcf2fasta -f $ref -p ${name}_${it}_ -P 1 ${name}_${it}.vcf
		mv ${name}_${it}_unknown* ${name}_${it}_cons.fasta

		### create vcf with lofreq
		rm ${name}_${it}_lofreq.vcf
		lofreq call -f $ref -o ${name}_${it}_lofreq.vcf ${name}_${it}.bam
	
		### calculate depth, run gap_cons.py
		samtools depth ${name}_${it}.bam > ${name}_${it}.depth
		$script_dir/gapcons.py ${name}_${it}_cons.fasta ${name}_${it}.depth $min_cov
	
		### remove temporary files of this iteration
		rm ${name}_${it}.sam
		rm ${name}_${it}.vcf
		rm smalt_index.*
		rm *.fasta.fai
		
		### new ref for next iteration
		ref=${name}_${it}_cons.fasta
		((it+=1))
	done
	
	### remove temporary files
	rm ${name}_reads.fastq
	
done

### Run spreadvcf.R
#wd=$(pwd)
#Rscript $script_dir/spreadvcf.R $wd