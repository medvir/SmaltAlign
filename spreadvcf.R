library(tidyverse)

path = "/Volumes/huber.michael/HCV/experiments/161229/"
files = list.files(path, pattern = "3_lofreq.vcf")

INFO_lofreq = c("DP", "AF", "SB", "DP4")

for (i in files) {
        name_i = gsub("_lofreq.vcf", "", i)
        vcf_file = paste0(path, name_i, "_lofreq.vcf")
        cons_file = paste0(path, name_i, "_cons_gaps.fasta")
        cov_file = paste0(path, name_i, "_cov.list")
        
        vcf_data =  read.table(vcf_file, quote="\"")
        colnames(vcf_data) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
        vcf_data = vcf_data %>%
                separate (INFO, INFO_lofreq, sep = ";", extra = "drop") %>%
                select(POS, REF, ALT, DP, AF)
        vcf_data$AF = gsub("AF=" ,"", vcf_data$AF)
        vcf_data$DP = gsub("DP=" ,"", vcf_data$DP)
    
        seq_data = data.frame(cons = unlist(strsplit(readLines(cons_file)[2], ""))) %>% mutate(POS = seq.int(nrow(.)))
        
        cov_data = data.frame(COV = readLines(cov_file)) %>% mutate(POS = seq.int(nrow(.)))
        
        comb_data = full_join(seq_data, vcf_data, "POS") %>%
                full_join(cov_data, "POS") %>%
                mutate(REF = as.character(REF)) %>%
                mutate(cons = as.character(cons)) %>%
                mutate(SEQ = ifelse(cons == REF, cons, paste0(REF, "*"))) %>%
                mutate(SEQ = ifelse(is.na(REF), cons, SEQ)) %>%
                #filter(POS >= 0 & POS <= 10000) %>% # filter for amplicon position
                select(POS, SEQ, ALT, AF, COV)
        
        DUPL = duplicated(comb_data$POS)
        comb_data = cbind(comb_data, DUPL)
        
        write.csv(comb_data, paste0(path, name_i, ".csv"))
}