# SmaltAlign
Initially, the script was used to make quick alignments of fastq reads against a reference using smalt, now it’s mainly used for HIV and HCV consensus generation for diagnostics.
`smaltalign.sh` is the main script. It does the following:
- sample reads with seqtk
- make de novo alignment of sampled reads with velvet
- align sampled reads and de novo contigs in triplicate against reference with smalt
- create consensus of iteration 1 with freebayes
- create vcf of iteration 1 with lofreq
- calculate depth of iteration 1 with samtools
- re-align sampled reads (without de novo contigs) against consensus of iteration 1
- create consensus, vcf and depth of iteration 2
- repeat for i iterations

All the necessary references are in the References directory.

### Installation of requirements.txt
To install all dependencies (listed above) you can use the command `conda install --file requirements.txt` (if [Bioconda](https://bioconda.github.io/) is installed).  
That's an easy way to assure a successful run of SmaltAlign.


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
`wts.R` is an R script to combine consensus sequence, variants and coverage for the last iteration of all `lofreq.vcf` files in a directory.
It saves a `WTS.fasta` file containing the consensus sequence with wobbles (at a certain threshold) and a `.csv` file  containing coverage and variant frequencies for every position.
The path to the working directory, the variant threshold and the minimal coverage have to be adapted manually in the first lines.
