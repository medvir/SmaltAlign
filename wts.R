### path to working directory
#path = commandArgs[1]
path = "/Volumes/huber.michael/Diagnostics/experiments/170309_HCV/"

### minority variant threshold (%)
variant_threshold = 5

### minimal coverage required (reads)
minimal_coverage = 3 

### libraries
library(tidyverse)
library(stringr)
library(seqinr)

files = list.files(path, pattern = "lofreq.vcf")
files = files[grep(paste0(max(str_sub(files, -12, -12)), "_lofreq.vcf"), files)] ### grep the N_lofreq.vcf files with the highest iteration

### functions
call_wobbles <- function(nuc_list) {
        ### function to call ambiguous nucleotides (wobbles)
        ### input is a list of nucleotides (nuc_list, length may be greater than 2)
        nuc_list = toupper(nuc_list) %>% unique() %>% .[!is.na(.)]  ### NAs in nuc_list are ignored (e.g. if position not covered, returns NA)
        if (length(nuc_list) == 0) {return(NA)}
        for (i in 1:length(nuc_list)) {
                if (nchar(nuc_list[i]) > 1) {return("insertion")}   ### returns "insertion" if nuc_list not composed of single characters
                if (!(is.element(nuc_list[i], c("A", "G", "C", "T", "N")))) {return("unkown_nucleotide")}   ### returns "unkown_nucleotide" if characters other than A, G, C, T, N in nuc_list
                }
        if ("N" %in% nuc_list | length(nuc_list) == 4) {return("N")}  ### returns "N" if at least one N or all four nucleotides included in nuc_list
        if (length(nuc_list) == 1) {return(nuc_list[1])}  ### one single nucleotide
        if (length(nuc_list) == 2) {  ### two different nucleotides
                IUPAC = data.frame(c("A", "R", "M", "W"), c("R", "G", "S", "K"), c("M", "S", "C", "Y"), c("W", "K", "Y", "T"))
                names(IUPAC) = c("A", "G", "C", "T")
                row.names(IUPAC) = c("A", "G", "C", "T")
                return(as.character(IUPAC[nuc_list[1], nuc_list[2]]))
        }
        if (length(nuc_list) == 3) {  ### three different nucleotides
                if (!("A" %in% nuc_list)) {return("B")}
                if (!("C" %in% nuc_list)) {return("D")}
                if (!("G" %in% nuc_list)) {return("H")}
                if (!("T" %in% nuc_list)) {return("V")}
                }
        }


### Loop over all files
for (i in files) {
        name_i = gsub("_lofreq.vcf", "", i)
        vcf_file = paste0(path, name_i, "_lofreq.vcf")
        cons_file = paste0(path, name_i, "_cons.fasta")
        depth_file = paste0(path, name_i, ".depth")
        
        if (class(try(read.table(vcf_file))) == "try-error") {next} ### next if vcf file is empty
        vcf_data = try(read.table(vcf_file, quote="\"")) %>%
                rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, QUAL = V6, FILTER = V7, INFO = V8) %>%
                separate (INFO, c("DP", "AF", "SB", "DP4"), sep = ";", extra = "drop") %>%
                select(POS, REF, ALT, DP, AF) %>%
                mutate(DP = gsub("DP=" ,"", DP)) %>%
                mutate(AF = round(as.numeric(gsub("AF=" ,"", AF))*100,1)) %>%
                filter(AF >= variant_threshold) ### filter for AF >= variant_threshold
    
        cons_data = data.frame(CONS = unlist(strsplit(readLines(cons_file)[-1], ""))) %>%
                mutate(POS = seq.int(nrow(.)))
            
        cov_data = read_delim(depth_file, "\t", col_names = FALSE, trim_ws = TRUE, col_types = "cii") %>%
                rename(POS = X2, COV = X3) %>%
                select(POS, COV)
        
        comb_data = full_join(cons_data, vcf_data, "POS") %>%
                full_join(cov_data, "POS") %>%
                #filter(POS <= 1000) %>% ### filter for positions
                mutate(REF = toupper(as.character(REF))) %>%
                mutate(ALT = toupper(as.character(ALT))) %>%
                mutate(CONS = toupper(as.character(CONS)))
        
        ### exit loop for this sample if alignment not correct
        if (!(all(comb_data$CONS == comb_data$REF, na.rm = TRUE))) {next} 
        
        ### call wobbles
        comb_data = comb_data %>%
                select(POS, CONS, ALT, AF, COV) %>%
                rename(REF = CONS) %>%
                mutate(WTS = apply(.[,c('REF', 'ALT')], 1, function(x) call_wobbles(c(x['REF'], x['ALT'])))) 
        
        ### resolve duplicate lines in vcf_data for two variants
        for (j in nrow(comb_data):2) {
                if (comb_data$POS[j] == comb_data$POS[j-1]) {
                        comb_data$WTS[j-1] = call_wobbles(c(comb_data$REF[j], comb_data$REF[j-1], comb_data$ALT[j], comb_data$ALT[j-1]))
                        comb_data$AF[j-1] = paste(comb_data$AF[j], comb_data$AF[j-1], sep = "/")
                        comb_data$ALT[j-1] = paste(comb_data$ALT[j], comb_data$ALT[j-1], sep = "/")
                        comb_data = comb_data[-j, ]
                }
        }
        
        ### apply minimal coverage threshold
        comb_data = comb_data %>%
                mutate(WTS = ifelse(is.na(COV), "-", ifelse(COV < minimal_coverage, "N", WTS)))
        
        ### write output files
        write.csv(comb_data, paste0(path, name_i,  "_", variant_threshold, ".csv"))
        write.fasta(paste(comb_data$WTS, collapse = ""), name_i, paste0(path, name_i, "_", variant_threshold, "_WTS.fasta"), open = "w", nbchar = 60)
}
