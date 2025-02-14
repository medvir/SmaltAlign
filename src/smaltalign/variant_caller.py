# -*- coding: utf-8 -*-
"""
The input is a raw fastq files directory and a reference sequence in a fasta file. 
For each fastq file in the input directory a new sample directory is created (based on the name of the input file) 
and all the results are stored in the corresponding directory.
By default the result directories are stored in the input dir unless the output_dir is given.
@status: development
"""
import sys
import os
import subprocess
import logging
import shutil
import glob

from smaltalign import construct_consensus

ITR = 4
velvet_kmer = 29
min_contig_lgth = 200

def run_child(cmd, exe='/bin/bash'):
    """Use subrocess.check_output to run an external program with arguments."""
    
    logging.info(cmd)
    try:
        output = subprocess.check_output(cmd, universal_newlines=True, shell=True, stderr=subprocess.STDOUT)
    
    except subprocess.CalledProcessError as ee:
        err_msg = 'error is %s \n' %(ee)
        logging.error(err_msg)
        output = None
    return output

def Quality_Control(input_file, path_name, min_qual_mean = 20):
    
    filtered_prefix = '%s_filtered' %path_name 
    filter_cmd = 'prinseq -fastq %s -min_qual_mean %d -log %s_prinseq.log \
    -out_good %s -out_bad %s_bad > %s_prinseq.err 2>&1' %(input_file, min_qual_mean, path_name, filtered_prefix, filtered_prefix, path_name)
    run_child(filter_cmd)
    #filtered_filename = filtered_prefix + '.fastq'
    return filtered_prefix
    
def main(input_path, reference_file, output_dir=None, default_reads=200000, 
         iterations=ITR, variant_th=15, min_cov=3, distant_ref=None, QC=False, min_qual_mean=20):
    
    # extract the base name
    dir_path = os.path.abspath(input_path)
    file_exts = ('*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz')
    all_input_files = []
    for f_ext in file_exts:
        file_path = os.path.join(dir_path, f_ext)
        all_input_files.extend(glob.glob(file_path))
    
    for input_file in all_input_files:

        ref_iterative_file = os.path.abspath(reference_file)
        
        if '_L001_R' in input_file:
            name = os.path.basename(input_file).split('_L001_R')[0]
        else:
            name = os.path.basename(input_file).split('.')[0]

        if output_dir == None:
            path_name = os.path.join(dir_path,name)
        else:
            output_abs_dir = os.path.abspath(output_dir)
            path_name = os.path.join(output_abs_dir,name)
        logging.info(path_name)
        
        if os.path.exists(path_name):
            shutil.rmtree(path_name)
        
        try:
            os.makedirs(path_name)
        except OSError as error:
            logging.error(error)
            continue
        
        if input_file.endswith('gz'):
            cmd_gz = 'gzip -d %s' %input_file
            run_child(cmd_gz)
            input_file = input_file.rstrip('.gz')
        
        reads_file = input_file
        path_name_prefix = os.path.join(path_name,name)
        if QC:
            logging.debug('filtering with prinseq')
            # can be added: -lc_method entropy -lc_threshold 70
            # -min_qual_mean 25
            filtered_prefix = Quality_Control(input_file, path_name_prefix, min_qual_mean)
            reads_file = '%s.fastq' %filtered_prefix
            
        # seqtk subsample the reads
        subsample_file = '%s_reads.fastq' %path_name_prefix
        if default_reads != -1:
            subsample_fastq_cmd = 'seqtk sample %s %i > %s ' %(reads_file, default_reads, subsample_file)
            run_child(subsample_fastq_cmd)
        else:
            shutil.copy(reads_file, '%s' %subsample_file )
            
        n_reads = 0 
        wc_fastq_cmd = 'wc -l %s' %(subsample_file)
        out = run_child(wc_fastq_cmd)
        n_reads = int(out.split()[0])//4
        
        # call velvet for assembly 
        velveth_cmd = 'velveth %s %d -fastq %s' %(path_name_prefix, velvet_kmer, subsample_file )
        run_child(velveth_cmd)
        
        velvetg_cmd = 'velvetg %s -min_contig_lgth %s' %(path_name_prefix, str(min_contig_lgth))
        run_child(velvetg_cmd)
        
        # convert contig.fasta file to fastQ file
        contig_path_name = os.path.join(path_name_prefix , 'contigs.fa')
        contig_fastq_path_name = os.path.join(path_name_prefix , 'contigs.fastq')
        convert_fa_fq_cmd = "seqtk seq -F '#' %s > %s" %(contig_path_name, contig_fastq_path_name)
        run_child(convert_fa_fq_cmd)
        
        # concat contig.fastq file to the reads file in triplicate
        cat_cmd = 'cat %s %s %s %s > %s_reads_contigs.fasta' %(subsample_file,  
                                                                contig_fastq_path_name, 
                                                                contig_fastq_path_name, 
                                                                contig_fastq_path_name, 
                                                                path_name_prefix)
        run_child(cat_cmd)
        
        for i in range(1,iterations+1):
            
            #sample $name, $n_sample reads, iteration $it
            logging.info('sample %s, %d reads, iteration %d \n'  %(name, n_reads, i ))
            
            
            # create iindex file for smalt-align
            smalt_index_cmd = 'smalt index -k 7 -s 2 %s_%d_smalt_index %s' %(path_name_prefix, 
                                                                             i, ref_iterative_file)
            run_child(smalt_index_cmd)
            
            samtools_index_cmd = 'samtools faidx %s' %(ref_iterative_file)
            run_child(samtools_index_cmd)
            
            # smalt align the reads to the reference
            smalt_align_cmd = ''
            if i == 1: # first iteration consider contigs from velvet 
                smalt_align_cmd = 'smalt map -n 24 -x -y 0.5 -f samsoft ' 
                smalt_align_cmd += '-o %s_%d.sam ' %(path_name_prefix, i)
                smalt_align_cmd += '%s_%d_smalt_index %s_reads_contigs.fasta' %(path_name_prefix, 
                                                                                i, path_name_prefix)
            else: 
                smalt_align_cmd = 'smalt map -n 24 -x -y 0.5 -f samsoft '
                smalt_align_cmd += '-o %s_%d.sam ' %(path_name_prefix, i)
                smalt_align_cmd += '%s_%d_smalt_index %s' %(path_name_prefix, 
                                                            i, subsample_file)
            
            run_child(smalt_align_cmd)
            
            # samtools to convert sam to bam and sort 
            samtools_sort_cmd = 'samtools view -Su %s_%d.sam | samtools sort -o %s_%d.bam' %(path_name_prefix, i, 
                                                                                             path_name_prefix, i)
            run_child(samtools_sort_cmd)
            samtools_index_cmd = 'samtools index %s_%d.bam'  %(path_name_prefix, i)
            run_child(samtools_index_cmd)
            
            # create consensus with freebayes
            freebayes_cmd = 'freebayes -f %s -p 1 %s_%d.bam > %s_%d.vcf' %(ref_iterative_file, 
                                                                           path_name_prefix, i, 
                                                                           path_name_prefix, i)
            run_child(freebayes_cmd)
        
            # check if vcf file has any mutations 
            muts_grep_cmd = 'grep -c -v "^#" %s_%d.vcf'  %(path_name_prefix, i)
            out = run_child(muts_grep_cmd)
            muts = int(out.split()[0])
            
            dst_file = '%s_%s_cons.fasta' %(path_name_prefix, i)
            
            # if no mutations 
            if muts == 0:
                logging.warnings('No variation found. \n')
                # copy the reference file to consensus file 
                shutil.copy(ref_iterative_file, dst_file)
            else:
                # use vcf2fasta to construct the consensus sequence
                vcf2fasta_cmd = 'vcf2fasta -f %s -p %s_%d_ -P 1 %s_%d.vcf' %(ref_iterative_file, 
                                                                             path_name_prefix, i, 
                                                                             path_name_prefix, i)
                run_child(vcf2fasta_cmd)
                outfile = '%s_%d_unknown*' %(path_name_prefix, i)
                src_path = glob.glob(outfile)[0]
                # rename the consensus sequence to ${name}_${it}_cons.fasta
                os.rename(src_path, dst_file) 
            
            # run lofreq to extract the vcf file 
            cpu_num = 1
            if os.cpu_count() > 1:
                cpu_num = os.cpu_count() // 2
            
            lofreq_vcf_file = '%s_%d_lofreq.vcf' %(path_name_prefix, i)
            lofreq_call_cmd = 'lofreq call-parallel --pp-threads %d ' %(cpu_num)
            lofreq_call_cmd += '-f %s -o %s ' %(ref_iterative_file, lofreq_vcf_file)
            lofreq_call_cmd += ' %s_%d.bam' %(path_name_prefix, i)
            run_child(lofreq_call_cmd)
            
            ## run commands to extract indels
            samtools_sort_bam_cmd = 'samtools sort %s_%d.bam -o %s_%d_sorted.bam'  %(path_name_prefix, i, 
                                                                                     path_name_prefix, i)
            run_child(samtools_sort_bam_cmd)
            
            lofreq_indelqual_file = 'lofreq indelqual --dindel -f %s %s_%d_sorted.bam  -o %s_%d_indel.bam' %(ref_iterative_file, 
                                                                                                             path_name_prefix, i, 
                                                                                                             path_name_prefix, i)
            run_child(lofreq_indelqual_file)
            
            samtools_ind_cmd = 'samtools index -b %s_%d_indel.bam' %(path_name_prefix, i)
            run_child(samtools_ind_cmd)
            
            lofreq_indel_file = 'lofreq call-parallel --pp-threads %d --call-indels -f %s -o %s_%d_lofreq_indel.vcf %s_%d_indel.bam' %(cpu_num, 
                                                                                                                                       ref_iterative_file, 
                                                                                                                                       path_name_prefix, i, 
                                                                                                                                       path_name_prefix, i)
            run_child(lofreq_indel_file)
            
            indels_grep_cmd = 'grep -c -v "^#" %s_%d_lofreq_indel.vcf'  %(path_name_prefix, i)
            out = run_child(indels_grep_cmd)
            indels = int(out.split()[0])
            if indels > 0:
                lofreq_filter_cmd = 'lofreq filter --only-indels -a 0.15 -v 10 --indelqual-thresh 20 -i %s_%d_lofreq_indel.vcf -o %s_%d_lofreq_indel_hq.vcf' %(path_name_prefix, i, 
                                                                                                                                                               path_name_prefix, i)
                run_child(lofreq_filter_cmd)
                logging.info('lofreq indel filter is done with %s' %lofreq_filter_cmd)
                
            #calculate depth for the covplot 
            depth_file = '%s_%d.depth' %(path_name_prefix, i)
            samtools_depth_cmd = 'samtools depth -d 1000000 %s_%d.bam ' %(path_name_prefix, i)
            samtools_depth_cmd += '> %s'  %(depth_file)
            run_child(samtools_depth_cmd)
            
            #clean_up()
            os.remove('%s_%d.sam' %(path_name_prefix, i))
            os.remove('%s_%d.vcf' %(path_name_prefix, i))
            index_fils = glob.glob('%s_%d_smalt_index.*' %(path_name_prefix, i))
            for f in index_fils:
                os.remove(f)
            
            cons_temp_file = glob.glob('%s_*_cons.fasta.fai' %(path_name_prefix))
            for f in cons_temp_file:
                os.remove(f)
            
            ref_iterative_file = dst_file
        
        os.remove('%s' %(subsample_file)) 
        os.remove('%s_reads_contigs.fasta' %(path_name_prefix))
        shutil.rmtree(path_name_prefix)
        
        #main(ref_file, lofreq_vcf_file, depth_file, output_file_str=None, VARIANT_TH=15, MINIMAL_COVERAGE=3, distant_ref=None):
        consensus_params_str = 'call consensus , reference file %s, variant_file %s, depth file %s, path_name_prefix %s, variant_th %s, minimum coverage %d, distant reference %s' %(
                                ref_iterative_file, lofreq_vcf_file, depth_file, path_name_prefix, variant_th, min_cov, distant_ref)
        #print(print_str)
        logging.info(consensus_params_str)
        construct_consensus.main(ref_iterative_file, lofreq_vcf_file, depth_file, 
                                 path_name_prefix, variant_th, min_cov, distant_ref)
        
        
if __name__ == "__main__":
    #reference_file = '../../../SARS_CoV_2_2.fasta'
    #input_dir = '../../../py_script_test'
    #output_dir = '../../../py_script_test/test_outdir'
    main(input_dir=sys.argv[1] , reference_file=sys.argv[2])
    
