# SmaltAlign

A consensus calling pipeline provided in Bash.

The pipeline can be used for variant calling analyses as well as to generate a consensus sequence based on the most common frequencies of the variant calling results.

`smaltalign.sh` does the following:
 1. Subsample reads with seqtk (optional).
 2. Make de novo alignment of sampled reads with velvet.
 3. If more than one reference is provided into `<reference_file>`, then it runs `select_ref_whole.py`, which selects the best reference sequence for each FASTQ file.
 4. Align sampled reads against reference with smalt. In the first iteration, to help "anchor" 
 5. Create consensus with freebayes.
 6. Create VCF (variant calling format) with lofreq.
 7. Calculate depth with samtools.
 8. If max number of iteration reached, call the final consensus sequence using final vcf file and the given ambiguity threshold. Otherwise, repeat from step 4.
 9. `cov_plot.R` plots the coverage of all iterations.
 10. `wts.R` combines consensus sequence, variants and coverage for the last iteration.

Multiple useful references are in the References directory. 
Further information about cov_plot.R and wts.R is provided below.

SmaltAlign also offers the option to proxess influenza sequences and analyse each of the segments. This is done with the `batch_influenza.sh` script as follows:
1. Iteration over all `.fastq.gz` files in the current directory.
2. Create a folder for each sample containing segment 1-8 subfolders.
3. If a reference folder containing a reference per segment is not provided or does not exist, then it runs `select_ref_segments.py`, which selects the best reference sequence for each segment from a Influenza reference database (selected sequences from the [NCBI Influenza Virus Database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database))
![IV-A references](References/genomes_query.png)
4. Using the best reference sequence, run `smaltalign.sh` for each segment, including `cov_plot.R` and `wts.R`.

## System requirements

#### Operating Systems

`smaltalign.sh` and `batch_influenza.sh` run on Linux. They have been tested on Linux Ubuntu 24.04.3 LTS.

#### Software Dependencies

The pipeline is designed to run within a Conda/Mamba environment.

The following specific versions (or higher) are recommended for stability as of 2026:

| Category | Tool | Version (Tested) | Purpose |
| :--- | :--- | :--- | :--- |
| **Search/alignment** | [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | 2.17.0 | `megablast` task for initial reference selection |
| **Mapping** | [SMALT](https://www.sanger.ac.uk/tool/smalt/) | 0.7.6 | Primary iterative genomic mapping |
| | [BWA](https://github.com/lh3/bwa) | 0.7.19 | Reference segment selection |
| **Assembly** | [Velvet](https://github.com/dzerbino/velvet) | 1.2.10 | *De novo* assembly for alignment seeding |
| **Variant Calling** | [LoFreq](https://csb5.github.io/lofreq/) | 2.1.5 | High-sensitivity variant and indel calling |
| | [FreeBayes](https://github.com/freebayes/freebayes) | 1.3.10 | Intermediate consensus generation |
| **Utilities** | [Samtools](http://www.htslib.org/) | 1.23 | BAM/SAM file manipulation and indexing |
| | [SeqKit](https://bioinf.shenwei.me/seqkit/) | 2.12.0 | FASTA/Q header cleaning and manipulation |
| | [Seqtk](https://github.com/lh3/seqtk) | 1.5 | Fastq subsampling |
| | [Bedtools](https://bedtools.readthedocs.io/) | 2.31.1 | Coverage statistics calculation |
| | [vcflib](https://github.com/vcflib/vcflib) | 1.0.3 | VCF to FASTA conversion utilities |
| | [bam-readcount](https://github.com/genome/bam-readcount) | 1.0.1 | Site-specific base frequency extraction |
| **Environment** | [Python](https://www.python.org/) | 3.12+ | Logic for reference selection and summaries |
| | [R](https://www.r-project.org/) | 4.3.3 | Coverage plotting and IUPAC consensus |

The pipeline relies on Python and R for data manipulation and visualization.

**Python Libraries:**
- **Pandas (>= 1.5.0):** Used for BLAST result analysis and tabular summaries.
- **Biopython (>= 1.81):** Used for all sequence I/O operations.

**R Libraries:**
- **Tidyverse (>= 2.0.0):** includes `readr`, `dplyr`, `tidyr`, `ggplot2`, `stringr`, needed for data processing and plotting.
- **cowplot (>= 1.1.1):** Theme and formatting for coverage plots.
- **seqinr (>= 4.2-30):** Writing consensus sequences.

#### Hardware requirements

* CPU: Multi-core processor (at least 4 cores recommended). The pipeline automatically utilizes half of the available cores for parallel tasks.

* RAM: 8GB minimum; 16GB recommended for large datasets or simultaneous 8-segment influenza processing.

* Storage: At least 10GB of free space are recommended for temporary intermediate files (stored in /tmp).

## Installation Guide

For a complete installation, please follow this 3 steps:

1. Install Conda: you need to have [Conda](https://docs.conda.io/en/latest/) installed.
2. Create the environment using the provided environment.yml.

To ensure you have all dependencies needed for SmaltAlign, you can use the `environment.yml` file.  
With the command `conda env create -f <path>/environment.yml` you will create a copy of the smaltalign environment.  
You enter the environment with the command `conda activate smaltalign` (and leave it with `conda deactivate`).  
For more information visit following link to [Managing environments](https://conda.io/docs/using/envs.html).

3. Clone/Copy the scripts into a local directory and ensure they are executable:

`git clone https://github.com/medvir/SmaltAlign.git`
`chmod +x *.sh *.py *.R`

#### Typical Install Time
On a normal desktop computer, the installation is expected to take about 10 minutes (mostly dependent on internet speed for downloading packages).

## Demo
You can run a demo on a single mock sample using the core script, `smaltalign.sh`, and also using the `batch_influenza.sh` script.

To run `smaltaign.sh`, you will need a sample FASTQ (available in the folder `Demo`) and a single reference FASTA (available in the same `Demo` folder).
Then, you can copy the following command:
`./SmaltAlign/smaltalign.sh -r ./demo/Reference_sequences/demo_ref.fasta -o ./demo/example_output_smaltalign/ ./demo/demo_file.fastq.gz`

To run `batch_influenza.sh`, you will need the same FASTQ file, and a few reference FASTA files, one per segment (also available in the `Demo` folder).
Then, you can copy the following command:
`./SmaltAlign/batch_influenza.sh -s ./demo/ -o ./demo/example_output_batch_influenza/ -r ./demo/Reference_sequences/`

#### Expected output
The following files will be generated in ./demo_results:

* [Sample]_majority_cons.fasta: The final consensus sequence.

* [Sample]_[Thresh]_WTS.fasta: Consensus sequence (with Wobbles).

* [Sample]_[Thresh].csv: A table containing coverage and variant frequencies per position (also based on the requested threshold, varients with frequencies lower than the threshold will not be displayed).

* coverage.pdf: A visual plot showing how coverage evolved over the iterations.

#### Expected run time

On a normal desktop, `smaltalign.sh` is expected to take about 2-3 minutes for a single sample with 200,000 reads.
For `batch_influenza.sh`, the time will be about 15-20 minutes, as it runs `smaltalign.sh` for each of the 8 segments.

## Instructions for Use

### `smaltalign.sh` usage
smaltalign.sh -r <reference_file> [options] <fastq_file/directory>  
Options:  
	[-h or --help]  
	[-n or --numreads]  
	[-i or --iterations]  
	[-t or --varthres]  
	[-c or --mincov]  
	[-o or --outdir]  
	[-d or --indels]  
	[-m or --mergecov]  

Run `smaltalign.sh -h` for detailed information of the options and the default parameters.

##### Reference file
If you would like to run the script using one specific reference sequence, you can specify this reference into `<reference_file>` (in fasta format).  
However, if you are not sure about the closest reference sequence, you can also specify several probable reference sequences in `<reference_file>` (only one file with all the sequences in fasta format). In this case, the script chooses the closest reference sequence from the set of given sequences and constructs the consensus sequence based on that.

##### Fastq file / directory
If you would like to analyse only one FASTQ file at a time, you can specify it into `<fastq_file/directory>`. In this case, you can specify either the folder where the file is located, or the file itself.
However, if you would like to anlayse multiple FASTQ files at a time, you can specify the folder to all the FASTQ files into `<fastq_file/directory>`. In this case, you need to have all your FASTQ files in the same folder.

##### Number of reads
By default, the tool subsamples the sequencing reads to 200000. However, this number can be changed by the users.  
If you would like to analyse all sequencing reads without subsampling, then use -n "all".

## `batch_influenza.sh` usage

batch_influenza.sh [options]  
Options:  
	[-h or --help]  
	[-r or --refdir]  
	[-s or --sampledir]  
	[-n or --numreads]  
	[-i or --iterations]  
	[-t or --varthres]  
	[-c or --mincov]  
	[-o or --outdir]  

Run `batch_influenza.sh -h` for detailed information of the outputs and the default parameters.

## wts.R
`wts.R` is an R script to combine consensus sequence, variants and coverage for the last iteration of all `lofreq.vcf` files in a directory.  
It saves a `_x_WTS.fasta` file containing the consensus sequence with wobbles (at a certain threshold x) and a `.csv` file  containing coverage and variant frequencies for every position.  
The variant threshold and the minimal coverage have to be provided by the user either when using `smaltalign.sh` or `batch_influenza.sh`.

## cov_plot.R
`cov_plot.R` is an R script to plot and save the coverage of all iterations of all `.depth` files in the working directory.

## Important Workflow Notes

* Iterative Refinement: The pipeline maps, calls variants, and updates the reference $N$ times (default is 4). This reduces "reference bias" and is crucial for highly variable viral samples.
* Temporary Files: The script uses /tmp for processing. If a run crashes, check /tmp/SmaltAlignOutputs_XXXXX to clear space.

# Contributions
- Judith Bergadà-Pijuan*
- Maryam Zaheri
- Stefan Schmutz
- Osvaldo Zagordi
- Michael Huber**

*maintainer ; **group leader
