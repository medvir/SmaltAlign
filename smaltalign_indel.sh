#!/bin/bash

### defaults
script_dir=$( dirname "$(readlink -f "$0")" )
n_reads=1000000
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
		ref="$2"
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

#if [[ $sample_dir =~ \.gz$ ]]; then
#    echo "counting number of reads"
#    lines=`zcat $sample_dir|wc -l`
#    max_seq=`expr $lines / 4`
#    n_reads=$max_seq
#elif [[ $sample_dir =~ \.fastq$ ]]; then
#    echo "counting number of reads"
#    lines=`cat $sample_dir|wc -l`
#    max_seq=`expr $lines / 4`
#    n_reads=$max_seq
#fi
#echo "max n_reads is $n_reads"

### print arguments
echo -e 'sample_dir: ' $sample_dir
echo -e 'script_dir: ' $script_dir
echo -e 'ref: ' $ref
echo -e 'n_reads: ' $n_reads
echo -e 'iterations: ' $iterations


### loop over list of files to analyse
if [ -d $sample_dir ]; then list=$(ls $sample_dir | grep .fastq); else list=$sample_dir; fi
for i in $list; do
	name=$(basename $i | sed 's/_L001_R.*//' | sed 's/.fastq.gz//'| sed 's/.fastq//')

	### sample reads with seqtk
	seqtk sample $i $n_reads > ${name}_reads.fastq
	n_sample=$(wc -l ${name}_reads.fastq | cut -f 1 -d " ")
	n_sample=$(($n_sample / 4))

	### de novo aligment
    #start=`date +%s`
	echo
	echo sample $name, $n_sample reads, de-novo alignment
	echo "*******************************************************************************"
	velveth ${name} 29 -fastq ${name}_reads.fastq
	velvetg ${name} -min_contig_lgth 200

	### fake fastq of contigs and cat to reads in triplicate
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${name}/contigs.fa | awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+"$1"\n"gensub(/./, "I", "g", $2)}' > ${name}/contigs.fastq
	cat ${name}_reads.fastq ${name}/contigs.fastq ${name}/contigs.fastq ${name}/contigs.fastq > ${name}_reads_contigs.fasta
    #end=`date +%s`
    #runtime=$((end-start))
    #echo "Assmebling runtime is $runtime"


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
		#start=`date +%s`
        smalt index -k 7 -s 2 ${name}_${it}_smalt_index $ref
		samtools faidx $ref

		if [ "$it" -eq 1 ]; then
			smalt map -n 24 -x -y 0.5 -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads_contigs.fasta
		else
			smalt map -n 24 -x -y 0.5 -f samsoft -o ${name}_${it}.sam ${name}_${it}_smalt_index ${name}_reads.fastq
		fi

		samtools view -Su ${name}_${it}.sam | samtools sort -o ${name}_${it}.bam
		samtools index ${name}_${it}.bam
        #end=`date +%s`
        #runtime2=$((end-start))
        #echo "Alignment smalt runtime is $runtime2"

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
		samtools sort ${name}_${it}.bam > ${name}_${it}_sorted.bam
		
		lofreq indelqual --dindel -f $ref ${name}_${it}_sorted.bam  -o ${name}_${it}_indel.bam
		samtools index -b ${name}_${it}_indel.bam
		lofreq call-parallel --pp-threads $cores --call-indels -f $ref -o ${name}_${it}_lofreq_indel.vcf ${name}_${it}_indel.bam

		if [ $(grep -v "^#" ${name}_${it}_lofreq_indel.vcf | wc -l) -gt 0 ] ; then
            lofreq filter --only-indels -a 0.15 -v 10 --indelqual-thresh 20 -i ${name}_${it}_lofreq_indel.vcf -o ${name}_${it}_lofreq_indel_temp_hq.vcf
        	echo lofreq filter is done ${name}_${it}_lofreq_indel_temp_hq.vcf 
        	if [ $(grep -v "^#" ${name}_${it}_lofreq_indel_temp_hq.vcf | wc -l) -gt 0 ] ; then
                cp ${name}_${it}_lofreq_indel_temp_hq.vcf ${name}_lofreq_indel_hq.vcf
                echo lofreq_indel_hq.vcf file is created. 
            fi
        fi
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