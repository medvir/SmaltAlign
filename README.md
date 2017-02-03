# SmaltAlign
Quick alignment of reads against a given reference using smalt. `smaltalign.sh` is the main script.

## Usage
	smaltalign.sh -r <reference_file> [options] <fastq_file/directory> 
	
	OPTIONS
	-r       reference
	-n INT   number of reads (default 50'000)
	-c INT   minimal coverage (default 3)
	-i INT   iterations (default 2)

### gapcons.py
`gapcons.py` is a Python script automatically called by smaltalign.sh to replace a position in a given sequence with "-" for coverage zero and with "N" for coverage less than n reads.
Arguments are 1) ref_file in fasta, 2) the samtools depth file, 3) minimal coverage n

### spreadvcf.R
`spreadvcf.R` is an R script to combine vcf and consenus for Cyril

### cov_plot.R
`cov_plot.R` is an R script to plot the coverage