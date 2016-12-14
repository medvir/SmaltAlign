#!/usr/local/bin/python3.4
'''
Replace reference sequence with "-" for zero coverage and with "N" for coverage less than n reads
Arguments are 1) ref_file in fasta, 2) the samtools depth file, 3) minimal coverage n.
'''

import sys
import os
import re

def cov_thres(seq_list, cov_list, n):

    ### assert that seq_list is longer than cov_list
    assert len(seq_list) >= len(cov_list)
        
    ### replace reference with "-" for coverage 0 and with "N" for coverage less than n reads
    for i, c in enumerate(cov_list):
        if c == 0:
            seq_list[i] = "-"
        
        elif c < n:
            seq_list[i] = "N"
    
    ### replace reference with "-" from end of coverage to end of reference
    for i in range(len(cov_list), len(seq_list)):
        seq_list[i] = "-"
    
    return(seq_list)


def depth_to_complete_list(depth_file):
    
    cov = {}
    for line in open(depth_file, 'r'):
        name, pos, v = line.strip().split('\t')
        cov[int(pos)] = int(v)
    
    cov_list = []
    for i in range(1, max(cov.keys()) + 1):
        cov_list.append(cov.get(i, 0))
    
    return(cov_list)        


def fasta_to_list(seq_file):

    ### write ref_file file into a list (seq_list)
    with open(seq_file, 'r') as file:
        seq = file.read().splitlines(True)[1:]
    seq = ''.join(seq)
    seq = seq.replace('\n', '')
    seq_list = list(seq)
    
    return(seq_list)


def main():
    
    seq_file = sys.argv[1]
    depth_file = sys.argv[2]
    n = int(sys.argv[3])
    out_file_gaps = seq_file.replace(".fasta", "_gaps.fasta")
    out_file_no_gaps = seq_file.replace(".fasta", "_no_gaps.fasta")
    out_file_cov_list = seq_file.replace("_cons.fasta", "_cov.list")
    name = seq_file.replace(".fasta", "")

    ### create list of the sequence
    seq_list = fasta_to_list(seq_file)
    
    ### create list with the complete coverage
    cov_list = depth_to_complete_list(depth_file)
    
    ### set coverage threshold
    seq_list = cov_thres(seq_list, cov_list, n)
        
    ### Join to string, remove "-", leave "N"
    seq = (''.join(seq_list))
    print(seq)
    
    ### write gap seq file into file
    with open(out_file_gaps, 'w') as file:
        file.write(">%s\n%s\n" %(name, seq))
        file.close()
        
    ### write gap seq file into file
    #seq = seq.replace("-", "")
    #with open(out_file_no_gaps, 'w') as file:
        #file.write(">%s\n%s\n" %(name, seq))
        #file.close()
		
	### write complete coverage file into file
    #with open(out_file_cov_list, 'w') as file:
        #for item in cov_list:
            #file.write("%s\n" % item)
		#file.close()


if __name__ == '__main__':
    main()