# SmaltAlign
Quick alignment of reads against a given reference using smalt. `smaltalign.sh` is the main script.

## Usage
	smaltalign.sh -r <reference_file> [options] <fastq_file/directory> 
	
	OPTIONS
	-r       reference
	-n INT   number of reads (default 50'000)
	-c INT   minimal coverage (default 4)
	-i INT   iterations (default 2)

### batch.sh
For running multiple samples in the current working directory with different references in one batch.

### gapcons.py
`gapcons.py` is a python script automatically called by smaltalign.sh to replace a position in a given sequence with "-" for coverage zero and with "N" for coverage less than n reads.
Arguments are 1) ref_file in fasta, 2) the samtools depth file, 3) minimal coverage n

### cov_plot.R
`cov_plot.R` is an R script to plot the coverage.
The path to the working directory has to be adapted in the first line.

### spreadvcf.R
`spreadvcf.R` is an R script to combine vcf and consenus sequence and make wobble bases for Cyril.
The path to the working directory has to be adapted in the first line.