# SmaltAlign
Quick alignment of reads against a given reference using smalt. `smaltalign.sh` is the main script.

## Usage
	smaltalign.sh -r <reference_file> [options] <fastq_file/directory> 
	
	OPTIONS
	-r       reference
	-n INT   number of reads (default 50'000)
	-i INT   iterations (default 2)

### batch.sh
For running multiple samples in the current working directory with different references in one batch.

### cov_plot.R
`cov_plot.R` is an R script to plot the coverage.
The path to the working directory has to be adapted in the first line.

### wts.R
`wts.R` is an R script to combine consenus sequence, variants and depth. It creates an output sequence with wobbles (at a certain threshold) 
The path to the working directory has to be adapted in the first line.