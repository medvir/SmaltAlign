# SmaltAlign
Initially used to make quick alignments of fastq reads against a reference using smalt, now mainly used for HIV and HCV consensus generation.
`smaltalign.sh` is the main script. It does the following:
- sample reads
- make de novo alignement 
- align reads and de novo contigs against reference
- create consensus and vcf round 1
- align reads (without de novo contigs) against consensus of round 1
- create consensus and vcf round 2
- ...
All the necessary references are in the References directory.

### Usage
	smaltalign.sh -r <reference_file> [options] <fastq_file/directory> 
	
	OPTIONS
	-r       reference
	-n INT   number of reads (default 200'000)
	-i INT   iterations (default 4)

### batch.sh
Used to run multiple samples in the current working directory with different references in one batch.

### cov_plot.R
`cov_plot.R` is an R script to plot and save the coverage of all iterations of all `.depth` files in the working directory.
The path to the working directory has to be adapted manually in the first line.

### wts.R
`wts.R` is an R script to combine consenus sequence, variants and coverage for the last iteration of all `lofreq.vcf` files in a directory.
It saves a `WTS.faste` file containging the consensus sequence with wobbles (at a certain threshold) and a `.csv` file  containing coverage and variant frequencies for every position.
The path to the working directory, the variant threshold and the minimal coverage have to be adapted manually in the first lines.