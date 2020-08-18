#!/opt/miniconda/bin/python3
"""Summarise SmaltAlign results from an influenza run.

Run on SmaltAlign results and extract base counts with bam-readcount. Then, extract coverage summary statistics
with bedtools.
"""
import sys
import os
import glob
import re
import shlex
import subprocess

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def mean_coverage(bamfile):
    """Run samtools depth to compute mean and stdev of coverage"""
    from math import sqrt
    cml = shlex.split('samtools depth -a -d 100000 %s' % bamfile)
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    with proc.stdout as handle:
        m = var = 0.0
        for l in handle:
            pos, c = (int(i) for i in l.strip().split()[1:])
            m += c
            var += c * c
    return round(m / pos, 2), round(sqrt(var / (pos - 1)), 2)

def coverage_info(bamfile, thresh=50):
    """Run bedtools to extract coverage information."""
    cml = shlex.split("bedtools genomecov -ibam %s -max %d" % (bamfile, thresh))
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    with proc.stdout as handle:
        for l in handle:
            if not l.startswith('genome'):
                continue
            lsp = l.strip().split()
            if int(lsp[1]) == thresh:
                return float(l.strip().split()[-1])
    return 0.0

def get_bases(reffile, bamfile, name=''):
    """Run bam-readcount, parse the results, and return a dictionary with base information."""
    info_dict = {
        'refname': [],
        'position': [],
        'refbase': [],
        'coverage': [],
        '=': [],
        'A': [],
        'C': [],
        'G': [],
        'T': [],
        'N': []
    }
    cml = shlex.split('bam-readcount -f %s -w 1 %s %s' % (reffile, bamfile, name))
    proc = subprocess.Popen(cml, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    with proc.stdout as handle:
        for i, line in enumerate(handle):
            chrom, pos, ref_b, depth = line.split('\t')[:4]
            info_dict['refname'].append(chrom)
            info_dict['position'].append(pos)
            info_dict['refbase'].append(ref_b)
            info_dict['coverage'].append(depth)
            for bi in line.split('\t')[4:]:
                base_info = bi.strip().split(':')
                try:
                    info_dict[base_info[0]].append(int(base_info[1]))
                except KeyError:
                    print('%s entities %s found at pos %s' % (base_info[1], base_info[0], pos), file=sys.stderr)
            if i == 4:
                pass
    return info_dict


def main():
    """What the main does."""
    rundir = sys.argv[1]
    run_name = os.path.split(rundir.rstrip('/'))[1]
    try:
        os.mkdir('tables-%s' % run_name)
    except OSError:
        sys.exit('Directory exists. Delete/rename old one.')
    molis_vec = []
    sample_vec = []
    segment_vec = []
    cov_20 = []
    cov_1k = []
    ave_cov = []
    stdev_cov = []
    # Run on all segments of all samples
    for segmentdir in glob.glob('%s/1000*_S*_L001_R1_001/segment*/' % rundir):
        print('## running %s ##' % segmentdir, file=sys.stderr)
        molis, sample, segment = re.search(r"\d*_[^/]*/(\d{10}).*(S\d*)_.*/segment(\d*)/", segmentdir).group(1, 2, 3)
        # extract iteration 4, skip if absent
        try:
            ref = glob.glob('%s/*4_cons.fasta' % segmentdir)[0]
            bam = glob.glob('%s/*4.bam' % segmentdir)[0]
        except IndexError:
            print('files not found in %s' % segmentdir, file=sys.stderr)
            continue
        refname = list(s.id for s in SeqIO.parse(ref, 'fasta'))[0]
        baminfo = get_bases(ref, bam, refname)
        # write base information to pandas dataframe, manipulate
        df = pd.DataFrame(baminfo)
        l = df.shape[0]
        df['molis'] = [molis] * l
        df['segment'] = [segment] * l
        df.drop(columns=['='], inplace=True)
        # extract nucleotides with max frequency to write consensus sequence
        nucleotides = df.loc[:, ['A', 'C', 'G', 'T', 'N']]
        nucleotides['MAX'] = nucleotides.idxmax(axis=1)
        seq = ''.join(nucleotides['MAX'].tolist())
        sr = SeqRecord(Seq(seq), id='%s-%s' % (molis, segment), description='')
        SeqIO.write(sr, 'tables-%s/%s-%s.fasta' % (run_name, molis, segment), 'fasta')
        # save table
        filename = 'tables-%s/%s-%s.csv' % (run_name, molis, segment)
        df = df.set_index('position')
        df.to_csv(filename)
        molis_vec.append(molis)
        sample_vec.append(sample)
        segment_vec.append(segment)
        cov_20.append(coverage_info(bam, 20))
        cov_1k.append(coverage_info(bam, 1000))
        m, s = mean_coverage(bam)
        ave_cov.append(m)
        stdev_cov.append(s)

        print('###############################', file=sys.stderr)
    # save coverage info to file
    cov_df = pd.DataFrame({'molis': molis_vec,
                           'segment': segment_vec,
                           'cov_20': cov_20,
                           'cov_1k': cov_1k,
                           'ave_cov': ave_cov,
                           'stdev_cov': stdev_cov})
    cov_df.to_csv('tables-%s/coverage_info.csv' % rundir, index=False)
if __name__ == '__main__':
    main()
