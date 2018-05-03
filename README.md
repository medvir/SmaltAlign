# SmaltAlign

[![Build Status](https://travis-ci.org/medvir/SmaltAlign.svg?branch=master)](https://travis-ci.org/medvir/SmaltAlign)

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

### Use conda environment from file
To ensure you have all dependencies needed for SmaltAlign installed you can use the `environment.yml` file.  
First you need to have [Conda](https://conda.io/docs/install/quick.html) installed).  
With the command `conda env create -f <path>/environment.yml` you will create a copy of the SmaltAlign environment.  
You enter the environment with the command `source activate SmaltAlign` (and leave it with `source deactivate`).  
For more information visit following link to [Managing environments](https://conda.io/docs/using/envs.html).


### Usage
	smaltalign.sh -r <reference_file> [options] <fastq_file/directory>

	OPTIONS
	-r       reference
	-n INT   number of reads (default 200'000)
	-i INT   iterations (default 4)

### batch.sh
Used to run multiple samples in the current working directory with different references in one batch.  
To analyse the results of a Diagnostic sequencing run following steps need to be done:
* create a new folder in `/data/Diagnostics/experiments/` with the date of the sequencing run (start-date, yymmdd)
* in that new folder create links to the .fastq files you want to analyse (`ln -sv`) and copy the `SampleSheet.csv` of that run
* copy the `batch.sh` file into that new folder
* add the filenames (you can use `sampleID_to_filename.xltx`) to the empty virus arrays in `batch.sh` separated by a new line (works if you copy from the excel file)
* activate SmaltAlign environment (`source activate SmaltAlign`)
* execute `./batch.sh`

After the iterations of every sample the .R scripts are executed and the progress is visible in `coverage.pdf`.

### cov_plot.R
`cov_plot.R` is an R script to plot and save the coverage of all iterations of all `.depth` files in the working directory.

### wts.R
`wts.R` is an R script to combine consensus sequence, variants and coverage for the last iteration of all `lofreq.vcf` files in a directory.
It saves a `_x_WTS.fasta` file containing the consensus sequence with wobbles (at a certain threshold x) and a `.csv` file  containing coverage and variant frequencies for every position.
The the variant threshold and the minimal coverage have to be adapted manually in the first lines.