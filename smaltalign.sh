#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
n_reads=50000
min_cov=3
iterations=3
readlength=150 ### can be used to trim reads

### arguments
if [ $# == 0 ]; then
	echo
	echo 'Usage: smaltalign.sh -r <reference_file> [options] <fastq_file/directory> '
	echo 'Options:'
	echo '  -r       reference'
	echo '  -n INT   number of reads (default 50000)'
	echo '  -c INT   minimal coverage (default 3)'
	echo '  -i INT   iterations (default 3)'
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
			gunzip -c $i | seqtk sample - $n_reads | /usr/local/seqtk/seqtk trimfq -L $readlength - > ${name}_reads.fastq
		else
			seqtk sample $i $n_reads | /usr/local/seqtk/seqtk trimfq -L $readlength - > ${name}_reads.fastq
		fi
		
	n_sample=$(wc -l ${name}_reads.fastq | cut -f 1 -d " ")
	n_sample=$(($n_sample / 4))
	
	### de novo aligment
	echo
	echo sample $name, $n_sample reads, de-novo alignment 
	echo "************************************************************"
	velveth ${name} 29 -fastq ${name}_reads.fastq
	velvetg ${name} -min_contig_lgth 100
	
	### fake fastq of contigs and cat to reads in triplicate
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${name}/contigs.fa | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+"$1"\n"gensub(/./, "I", "g", $2)}' > ${name}/contigs.fastq
	cat ${name}_reads.fastq ${name}/contigs.fastq ${name}/contigs.fastq ${name}/contigs.fastq > ${name}_reads_contigs.fasta
		
	it=1
	while [ "$it" -le "$iterations" ]; do		
		echo
		echo sample $name, $n_sample reads, iteration $it
		echo "************************************************************"
				
		### align with smalt to reference, reads either ${name}_reads.fastq or ${name}_reads_contigs.fasta
		smalt index -k 7 -s 2 ${name}_${it}_smalt_index $ref
		if [ "$it" -eq 1 ]; then
			smalt map -n 28 -x -y 0.5 -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads_contigs.fasta
		else
			smalt map -n 28 -x -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads.fastq
		fi
		
		samtools view -Su ${name}_${it}.sam | samtools sort - ${name}_${it}
		samtools index ${name}_${it}.bam	

		### create consensus with freebayes
		freebayes -f $ref -p 1 ${name}_${it}.bam > ${name}_${it}.vcf	
		vcf2fasta -f $ref -p ${name}_${it}_ -P 1 ${name}_${it}.vcf
		mv ${name}_${it}_unknown* ${name}_${it}_cons.fasta

		### create vcf with lofreq
		rm -f ${name}_${it}_lofreq.vcf
		lofreq call -f $ref -o ${name}_${it}_lofreq.vcf ${name}_${it}.bam
	
		### calculate depth, run gap_cons.py
		samtools depth ${name}_${it}.bam > ${name}_${it}.depth
		$script_dir/gapcons.py ${name}_${it}_cons.fasta ${name}_${it}.depth $min_cov
	
		### remove temporary files of this iteration
		rm ${name}_${it}.sam
		rm ${name}_${it}.vcf
		rm ${name}_${it}.depth
		rm ${name}_${it}_smalt_index.*
		rm ${name}_*_cons.fasta.fai
		
		### new ref for next iteration
		ref=${name}_${it}_cons.fasta
		((it+=1))
	done
	
	### remove temporary files
	#rm ${name}_reads.fastq
	rm ${name}_reads_contigs.fasta
	rm -rf ${name}
done