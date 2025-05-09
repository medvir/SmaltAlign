#!/opt/miniconda/bin/python3
import sys
import subprocess

from textwrap import fill
import pandas as pd
from Bio import SeqIO

readfile = sys.argv[1]
outdir = sys.argv[2]
refdir = sys.argv[3]

allrefs = dict([(s.id.split('_')[0], str(s.seq))
                for s in SeqIO.parse('%s/References/flugenomes_nonmixed.fasta' % refdir, 'fasta')])

# index flugenomes.fasta
cml = 'bwa index %s/References/flugenomes_nonmixed.fasta' % refdir
subprocess.call(cml, shell=True)
# align against all genomes
cml = 'bwa mem -t 24 %s/References/flugenomes_nonmixed.fasta %s | samtools view -F 4 > %s/aln.sam' % (refdir, readfile, outdir)
subprocess.call(cml, shell=True)
# extract accession number, segment, serotype
cml = 'cut -f 3 %s/aln.sam | cut -d "_" -f 1-3 | tr -d ">" | tr "_" "\t" > %s/ref.tsv' % (outdir, outdir)
subprocess.call(cml, shell=True)

# manipulate with pandas to find, for each segment, the sequence with most hits
df = pd.read_table('%s/ref.tsv' % outdir, names=['accn', 'segment', 'serotype'])
count_ref = df.groupby(['segment', 'accn', 'serotype']).size()
c = count_ref.reset_index(name='counts').sort_values(['segment', 'counts'], ascending=[True, False])
c.to_csv('%s/counts.tsv' % outdir, index=False, sep='\t')
print(c.groupby('segment').head(3))

for segment in range(1, 9):
    counts = c[c['segment'] == segment]
    if counts.counts.sum() < 200 or counts.empty:
        print(segment, 'not enough')
        continue
    best_acc = counts.accn.tolist()[0]
    print(segment, best_acc)
    best_seq = allrefs[best_acc]
    with open('%s/segment-%d.fasta' % (outdir, segment), 'w') as h:
        h.write('>segment-%d-%s\n' % (segment, best_acc))
        h.write(fill(best_seq, width=80))
