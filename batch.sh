#!/bin/bash
	
n=50000 # number of reads (default 50000)
i=4     # iterations (default 4)

G1a=("1000310031-G1a-226-1_S11_L001_R1_001.fastq.gz
1000310742-G1a-226-2_S12_L001_R1_001.fastq.gz
1000311130-G1a-226-3_S13_L001_R1_001.fastq.gz
1000311660-G1a-226-4_S14_L001_R1_001.fastq.gz
1000312141-G1a-226-5_S15_L001_R1_001.fastq.gz
1000312484-G1a-226-6_S16_L001_R1_001.fastq.gz
1000338688-G1a-225-3_S1_L001_R1_001.fastq.gz
1000340531-G1a-229-6_S18_L001_R1_001.fastq.gz
1000340531-G1a-229-6_S4_L001_R1_001.fastq.gz
1000340531-G1a-229-N2-2-2-4_S8_L001_R1_001.fastq.gz")

G1l=("1000339557-G1l-225-7_S2_L001_R1_001.fastq.gz")

G1b=("1000340833-G1b-229-3_S5_L001_R1_001.fastq.gz
1000340833-G1b-229-7_S19_L001_R1_001.fastq.gz
1000340833-G1b-229-N2-3-1-4_S9_L001_R1_001.fastq.gz")

G3a=("")

G4d=("1000339981-G4d-229-1_S3_L001_R1_001.fastq.gz
1000339981-G4d-229-5_S17_L001_R1_001.fastq.gz
1000339981-G4d-229-N2-1-1-4_S7_L001_R1_001.fastq.gz
1000340953-G4d-229-4_S6_L001_R1_001.fastq.gz
1000340953-G4d-229-8_S20_L001_R1_001.fastq.gz
1000340953-G4d-229-N2-4-1-4_S10_L001_R1_001.fastq.gz")

HIV2=("")




script_dir=$( dirname "$(readlink -f "$0")" )

for s in $G1a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_1.fasta -n $n -i $i $s
done

for s in $G1b; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_1b_M1LE.fasta -n $n -i $i $s
done

for s in $G1l; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_1l_KC248196.fasta -n $n -i $i $s
done

for s in $G3a; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_3a_CB.fasta -n $n -i $i $s
done

for s in $G4d; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HCV_4d.fasta -n $n -i $i $s
done

for s in $HIV2; do
	${script_dir}/smaltalign.sh -r ${script_dir}/References/HIV-2A_ROD.fasta -n $n -i $i $s
done

#Rscript ${script_dir}/cov_plot.R ./
#Rscript ${script_dir}/wts.R ./