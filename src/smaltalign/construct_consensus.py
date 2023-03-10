#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@status: development
"""

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

from Bio.SeqIO import FastaIO
import os
import glob
import re
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from io import StringIO
import logging

d2a = {'AG': 'R', 'CT': 'Y', 'AC': 'M', 'GT': 'K', 'CG': 'S', 'AT': 'W',
       'ACT': 'H', 'CGT': 'B', 'ACG': 'V', 'AGT': 'D', 'ACGT': 'N'}

NC = {'A','C','G','T','N'}

wobbles = {v: k for k, v in d2a.items()}


def main(ref_file, lofreq_vcf_file, depth_file, output_file_str=None, VARIANT_TH=15, MINIMAL_COVERAGE=3, distant_ref=None):

    # read vcf file (from smaltalign lofreq output) 
    if not output_file_str:
        output_file_str = lofreq_vcf_file.split('_lofreq')[0] #lofreq_vcf_file.split('%d_lofreq' %max_itr)[0]
    
    logging.info('Output directory: %s' %output_file_str)
    
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    try:
        lofreq_vcf_pd = pd.read_csv(lofreq_vcf_file, sep='\t', comment='#', header=None,  names=columns)
    except FileNotFoundError:
        logging.error( '%s file not found' %lofreq_vcf_file)
        
    lofreq_vcf_pd[["DP", "AF", "SB", "DP4"]] = lofreq_vcf_pd['INFO'].str.split(';', 3, expand=True)
    lofreq_vcf_pd['AF'] = lofreq_vcf_pd['AF'].map(lambda x: round(float(x.lstrip('AF=')) * 100, 1))
    lofreq_vcf_pd['DP'] = lofreq_vcf_pd['DP'].map(lambda x: int(x.lstrip('DP=')))
    lofreq_vcf_pd['SB'] = lofreq_vcf_pd['SB'].map(lambda x: int(x.lstrip('SB=')))
    
    # filter for low allele frequencies
    vcf_pd = lofreq_vcf_pd[['POS', 'REF', 'ALT', 'DP', 'AF']].copy()
    if VARIANT_TH > 0: # not majority consensus
        vcf_pd = vcf_pd[vcf_pd['AF'] >= VARIANT_TH]
    
    vcf_pd['AF'] = vcf_pd['AF'].astype(str)
    
    aggregation_functions = {'REF': 'first', 'ALT': ','.join, 'DP': 'mean', 'AF': ','.join}
    vcf_pos_grouped_pd = vcf_pd.groupby(vcf_pd['POS']).aggregate(aggregation_functions).reset_index()
    
    # read the closest consensus sequence  (smaltalign freebayes + vcf2fasta output from the last iteration) 
    ref_seq = list(SeqIO.parse(ref_file, 'fasta'))[0]
    ref_pos_ls = [i+1 for i in range(len(ref_seq))]
    cons_pd = pd.DataFrame({'POS': ref_pos_ls, 'ref_cons': list(ref_seq)})

    # merge the consensus nucleotide, positions with vcf information
    cons_vcf_pd = pd.merge(cons_pd, vcf_pos_grouped_pd, on='POS', how="outer")
    cons_vcf_pd['REF'] = cons_vcf_pd[cons_vcf_pd['REF'].notnull()]['REF'].str.upper()  #['REF', 'ALT', 'ref_cons']
    cons_vcf_pd['ALT'] = cons_vcf_pd[cons_vcf_pd['ALT'].notnull()]['ALT'].astype(str).str.upper()
    cons_vcf_pd['ref_cons'] = cons_vcf_pd[cons_vcf_pd['ref_cons'].notnull()]['ref_cons'].str.upper()

    # read depth file (from smaltalign samtools output) 
    columns = ['ref','POS','COV']
    try:
        depth_pd = pd.read_csv(depth_file, sep='\t', comment='#', header=None,  names=columns)
    except FileNotFoundError:
        logging.error( '%s file not found' %depth_file)
        
    depth_pd.drop(['ref'] , axis=1, inplace=True)
    # merge the depth information with the consensus nucleotide, positions and vcf information
    cons_vcf_depth_pd = pd.merge(cons_vcf_pd, depth_pd, on='POS', how="outer")

    wobbles_dict = {}
    # for each position in the close reference sequence
    for index, row in cons_vcf_depth_pd.iterrows():
        nuc_ls = []
        
        # Notice here REF is in a seperate col and ALT in the format A,C in the different columns 
        cur_pos = row['POS']
        
        if pd.isna(row['REF']) and pd.isna(row['ref_cons']):
            #print('Error, ref and cons are both None.')
            logging.error("ref and cons are both None at position %d", cur_pos)
            continue
        
        if pd.isna(row['COV']):
            wobbles_dict[cur_pos] = '-'
            continue
        elif row['COV'] < MINIMAL_COVERAGE or 'N' in nuc_ls or len(nuc_ls) == 4:
       
            wobbles_dict[cur_pos] = 'N'
            continue
        
        if pd.isna(row['REF']):
            nuc_ls.append(row['ref_cons'])
        else:
            nuc_ls.append(row['REF'])
        
        if not pd.isna(row['ALT']):
            ALT_ls = row['ALT'].split(',')
            nuc_ls += ALT_ls 
            AF_ls = [float(i) for i in str(row['AF']).split(',')]
            assert len(ALT_ls) == len(AF_ls), logging.error("number of alleles are different from number of allele frequencies in position %d."  %(cur_pos))
            
            NCs_AF_dict = {}
            for i in range(len(ALT_ls)): 
                if ALT_ls[i] not in NCs_AF_dict:
                    NCs_AF_dict[ALT_ls[i]] =  AF_ls[i]
                    
            assert sum(AF_ls) <= 100, logging.error("sum(AF_ls) %d in position %d."  %(sum(AF_ls),cur_pos)) # if this happen print position
            
            ref_freq = 100 - sum(AF_ls) 
            NCs_AF_dict[nuc_ls[0]] = ref_freq
        
        if len(nuc_ls) == 0:
            
            logging.error('Error, No NC in position %d.' %(cur_pos))
            wobbles_dict[cur_pos] = '!'
            
        for n in nuc_ls:
            if len(n) > 1:
                logging.warning('unexpected insertion in position %d is %s.' %(cur_pos, n))
                #sys.exit()
                continue
            if n not in NC:
                logging.error('Error, nucleotide in position %d is %s' %(cur_pos, n))
                continue 
        
        if len(nuc_ls) == 1: # There is no allele at this position 
            wobbles_dict[cur_pos] = nuc_ls[0]
        elif len(nuc_ls) > 1: # there are allele at this position
            if ''.join(sorted(nuc_ls)) in d2a:
                if VARIANT_TH > 0: 
                    wobbles_dict[cur_pos] = d2a[''.join(sorted(nuc_ls))]
                else:
                    wobbles_dict[cur_pos] = max(NCs_AF_dict, key=NCs_AF_dict.get)    
            else: 
                logging.error('Error, nucleotides %s in position %d not in the d2a dictionary.' %(''.join(nuc_ls), cur_pos))
                sys.exit()
    
    cons_vcf_depth_pd['WTS'] = cons_vcf_depth_pd['POS'].map(wobbles_dict)
    cons_vcf_depth_pd.drop(['DP', 'REF'] , axis=1, inplace=True)
    
    o_filename = output_file_str + '_th%d.csv'%VARIANT_TH
    cons_vcf_depth_pd.to_csv(o_filename, index=False)
    
    amb_cons_seq = ''.join(cons_vcf_depth_pd["WTS"].fillna(''))
    amb_cons_id = output_file_str.split('/')[-1] + '_variant_th_' + str(VARIANT_TH)
    amb_record = SeqRecord(
        Seq(amb_cons_seq),
        id=amb_cons_id,
        description='',
    )
    
    final_cons_file_name = output_file_str + '_th%d_WTS.fasta'%VARIANT_TH
    output_handle = open(final_cons_file_name, "w")
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_record(amb_record)
    output_handle.close()
    
    if distant_ref:
        try:
            distant_ref_record = SeqIO.read(distant_ref, "fasta")
        except:
            sys.exit('Error reading distant reference file %s' %distant_ref)
        
        # Identical characters are given 2 points, 3 point is deducted for each non-identical character.
        # 5 points are deducted when opening a gap, and 2 points are deducted when extending it.
        # The parameters are taken from NCBI nucleotide BLAST 
        needle_cline = NeedleCommandline(auto=True, asequence = distant_ref, bsequence = final_cons_file_name, 
                                         gapopen=10, gapextend=0.5, stdout=True)
        stdout, stderr = needle_cline()
        
        alignments = list(AlignIO.parse(StringIO(stdout), "emboss"))
        
        distant_ref_record_aligned = alignments[0][0].seq # alignments[0].seqA #HXB2_record
        cons_aligned = alignments[0][1].seq #alignments[0].seqB #cons_record
        
        align_allpos_ls = []
        distant_ref_pos = 0
        cons_pos = 0
        cons_nc_pos = 0 
        distant_ref_nc_pos= 0 
        for i in range(len(distant_ref_record_aligned)):
            if cons_aligned[i] == '-':
                cons_pos = 0
            else:
                cons_nc_pos += 1
                cons_pos = cons_nc_pos
            
            if distant_ref_record_aligned[i] == '-':
                distant_ref_pos = 0
            else:
                distant_ref_nc_pos += 1
                distant_ref_pos = distant_ref_nc_pos
                
            align_allpos = [i, cons_pos, cons_aligned[i], distant_ref_pos, distant_ref_record_aligned[i]]
            align_allpos_ls.append(align_allpos)
        
        align_allpositions_df = pd.DataFrame(align_allpos_ls, columns = ['alignment_pos', 'cons_pos','cons', 'distant_ref_pos',distant_ref_record.id])
        
        o_filename = output_file_str + '_th%d_aligned'%VARIANT_TH + '.csv'
        align_allpositions_df.to_csv(o_filename, index=False)
        

if __name__ == "__main__":
    
    main(ref_file=sys.argv[1], lofreq_vcf_file=sys.argv[2], depth_file=sys.argv[3])
