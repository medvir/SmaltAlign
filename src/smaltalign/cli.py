#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@status: development
@author: maryamzaheri

"""

import argparse
import sys
from pkg_resources import (get_distribution, DistributionNotFound)

try:
    __version__ = get_distribution('smaltalign').version
except DistributionNotFound:
    __version__ = 'unknown'

# Parse command line
parser = argparse.ArgumentParser(prog="smaltalign", 
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Call consensus sequence using a raw fastq file and the given reference sequence 
in a fasta file and create a consensus sequence and a consensus information csv file. 
The input is a raw fastq files directory and a reference sequence in a fasta file. 
For each fastq file in the input directory a new sample directory is created (based on the name of the input file) 
and all the results are stored in the corresponding directory.
By default the result directories are stored in the input dir unless the output_dir is given.""",

epilog="""
example:
smaltalign -d '../../../py_script_test' -r references/SARS_CoV_2_2.fasta 

""")

#dest="f",
parser.add_argument("--fastq_dir", "-d", type=str, required=True,
    help="directory that has raw fastq files, A unique identifier from each fastq\n\
filename will be selected and the directory and all\n\
files generated from this file will have this identifier.\n\
If filename include 'XXX_L001_R*' then XXX will be the unique identifier\n\
else the name of the file before first '.' will be the unique identifier")

#metavar='REFERENCE'
parser.add_argument('--ref', '-r', type=str, required=True,
                    help='A close reference sequence or a close consensus\n\
sequence obtained iteratively using smaltalign.')
#metavar='PREFIX',
parser.add_argument('--output_dir', '-o',  type=str, default=None, 
                    help='The directory that will have all the results.\n\
If does not exist it will be created, if not given the input directory is used.')

#metavar='READS_NUMBER',
parser.add_argument('--reads_number', '-n', type=int, default=200000,
                    help='Subsample sequences. For no subsampling set this parameter to -1')

#metavar='ITERATIONS', 
parser.add_argument('--iterations', '-i', type=int, default=4 ,
                    help='Number of iterations the variants are called to construct the consensus sequence.\n\
For heighly mutated viruses it can be up to 4 and for viruses with\n\
low mutation rates such as SARS-CoV-2 it can be 2.')

#metavar='THRESHOLD',
parser.add_argument('--threshold', '-t',  type=float, default=15, 
                    help='variant threshold to construct consensus sequence, default is 15%%.')

#metavar='MINCOV',
parser.add_argument('--min_coverage', '-c',  type=int, default=3, 
                    help='minimum coverage to call variants in the final iteration.')

#metavar='DISTREFERENCE',
parser.add_argument('--distant_ref', '-f',  type=str, default=None, 
                    help='A distant reference sequence if positioning based on it is needed in the output file.')

parser.add_argument('-v', '--version', action='version',
                    version=__version__)

# Exit if no input file is provided
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()


def main(args=None):

    import logging
    import logging.handlers

    args = parser.parse_args()

    log_format = '%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s'
    logging.basicConfig(filename='smaltalign.log', level=logging.INFO, filemode='a',
                        format=log_format, datefmt='%Y/%m/%d %H:%M:%S')
    
    logging.info(' '.join(sys.argv))

    from smaltalign import variant_caller
    variant_caller.main(args.fastq_dir, args.ref, args.output_dir, 
                        args.reads_number, args.iterations, args.threshold, 
                        args.min_coverage, args.distant_ref)