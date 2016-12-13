# SmaltAlign
Quick alignment of reads agains a given reference using smalt.

## smaltalign.sh
`smaltalign.sh` is the script for HIV-1 env plasmid sequencing analysis.
Input: Fastq file or directory containing fastq files

### gapcons.py
`gapcons.py` is a script to replace a position in a given sequence with "-" for zero coverage and with "N" for coverage less than n reads.
Arguments are 1) ref_file in fasta, 2) the samtools depth file, 3) minimal coverage n

### Usage
	usage: smaltalign.sh -r <reference_file> [options] <fastq_file/directory> 
	
	OPTIONS
	-r       reference
	-n INT   number of reads (default 100'000)
	-c INT   minimal coverage (default 3)