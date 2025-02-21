#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:19:24 2024

This script select the most probable reference sequence from a FASTA file given 
a FASTQ file by performing the following steps:
1. Subsamples the FASTQ file to a specified number of reads.
2. Converts the subsampled reads to FASTA format.
3. Aligns the subsampled FASTA sequences to a reference FASTA file using BLAST.
4. Analyzes the BLAST results to identify the best matching sequence.
5. Extracts the best matching sequence from the reference FASTA file and 
writes it to a new FASTA file.

@author: maryamzaheri
"""
import argparse
import pandas as pd

from Bio import SeqIO

import subprocess
import sys
import logging
import logging.handlers
import shlex

def extract_fasta_by_id(input_fasta, sequence_id, output_fasta):
    """Extracts a FASTA sequence by ID and writes it to a new file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        sequence_id (str): ID of the sequence to extract.
        output_fasta (str): Path to the output FASTA file.
    """
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == sequence_id:
                with open(output_fasta, "w") as outfile:
                    SeqIO.write(record, outfile, "fasta")
                logging.info(f"Sequence {sequence_id} has been extracted and written to {output_fasta}.")
                return
        logging.info(f"Sequence ID {sequence_id} not found in the input file.")
        
def main(fastq_file, reference_fasta, subsample_num, outdir):
    """Performs a BLAST analysis and selects a best reference sequence.
    Args:
        fastq_file (str): Path to the FASTQ file.
        reference_fasta (str): Path to the reference FASTA file.
        subsample_num (int): Number of reads to subsample.
    """
    logging.info(f" fastq file: {fastq_file}, reference file: {reference_fasta}, subsampled sequences {subsample_num}")
    # subsample fastq file
    subsample_cmd = shlex.split('seqtk sample %s %d' % (fastq_file, subsample_num))
    subsample_process = subprocess.Popen(subsample_cmd, stdout=subprocess.PIPE)
    convert_cmd = shlex.split('seqtk seq -A -')
    with open('%s/subsampled_sequences.fasta' % outdir, 'w') as fasta_output:
        convert_process = subprocess.Popen(convert_cmd, stdin=subsample_process.stdout, stdout=fasta_output)
        subsample_process.stdout.close()
        convert_process.communicate()
    fasta_output.close()
    
    blast_output_file = '%s/blast_results.tsv' % outdir
    blast_cmd = 'blastn -task megablast -query %s/subsampled_sequences.fasta -outfmt 6' % outdir
    blast_cmd += ' -subject {} -out {}'.format(shlex.quote(reference_fasta),
                                                  blast_output_file)
    blast_process = subprocess.Popen(shlex.split(blast_cmd), universal_newlines=True)
    blast_process.wait()
    
    # read blast results and assigns to best subjects
    blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    blast_results = pd.read_csv(blast_output_file, names=blast_columns, delimiter="\t")
    if blast_results.empty:
        logging.warning('No read from the suggested virus was found')
        
    total_hits = blast_results.shape[0]
    unique_queries = len(set(blast_results['qseqid']))
    logging.info('BLAST queries: %d\tHits: %d', unique_queries, total_hits)
    
    hit_frequencies = {sequence_id: 0.0 for sequence_id in set(blast_results['sseqid'])}
    reference_frequencies = {sequence_id: 0.0 for sequence_id in set(blast_results['sseqid'])}
    
    grouped_results = blast_results.groupby(['qseqid'])
    for query_id, group in grouped_results:
        best_match = group[group['pident'] == group['pident'].max()]['sseqid']
        for match in best_match:
            hit_frequencies[match] += 1.0 / (unique_queries * len(best_match))
            reference_frequencies[match] += 1.0 / len(best_match)
    
    max_frequency = max(hit_frequencies.values())
    best_reference_id = max(reference_frequencies, key=reference_frequencies.get)
    logging.info(f"Best Reference ID {best_reference_id} with max frequency {max_frequency}.")
    
    extract_fasta_by_id(reference_fasta, best_reference_id, output_fasta = '%s/chosen_reference.fasta' % outdir)
    sorted_hit_frequencies = sorted(hit_frequencies.items(), key=lambda x:x[1], reverse=True)
    
    with open('%s/reference_freq.csv' % outdir, 'w') as ref_freq_output:
        for match, f in sorted_hit_frequencies:
            #if sorted_hit_frequencies[match]:
            print('%s,%5.4f' % (match, f), file=ref_freq_output)
        
    return best_reference_id, max_frequency 



if __name__ == "__main__":
    
    # parse command line
    parser = argparse.ArgumentParser()
    # First define all option groups
    group1 = parser.add_argument_group('Input file', 'Required input')
    group1.add_argument("-f", "--fastq", default="", type=str, dest="f",
                        help="input reads in fastq format")
    group1.add_argument("-r", "--fasta", default="", 
                        type=str, dest="r",
                        help="reference sequences in fasta format")
    group1.add_argument("-s", "--subsample", default="1000", type=int, dest="s",
                        help="number of subsampled reads")
    group1.add_argument("-o", "--outdir", default="./", type=str, dest="o",
                        help="path to the output directory")

    # exit so that log file is not written
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'
    logging.basicConfig(filename='%s/select_reference.log' % args.o, level=logging.INFO, format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    chosen_sequence_id, max_frequency = main(args.f, args.r, args.s, args.o)
