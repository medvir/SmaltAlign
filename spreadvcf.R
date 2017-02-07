path = "/Volumes/huber.michael/HCV/experiments/170202/"
threshold = 15 # minority variant threshold (%)
coverage = 3 # minimal coverage required (reads)

library(tidyverse)
library(stringr)
library(seqinr)

files = list.files(path, pattern = "lofreq.vcf")
files = files[grep(paste0(max(str_sub(files, -12, -12)), "_lofreq.vcf"), files)] # grep the lofreq files with the highest iteration
INFO_lofreq = c("DP", "AF", "SB", "DP4")

wobble <- function(l) {  ### input is a list of nucleotides
        l = toupper(l) %>% unique()
        if (!(all(l %in% c("A", "G", "C", "T", "N")))) {return("-")}  ### returns '-' if other characters than nucleotides
        if ("N" %in% l | length(l) == 4) {return("N")}
        if (length(l) == 1) {return(l[1])}
        if (length(l) == 2) {
                IUPAC = data.frame(c("A", "R", "M", "W"), c("R", "G", "S", "K"), c("M", "S", "C", "Y"), c("W", "K", "Y", "T"))
                names(IUPAC) = c("A", "G", "C", "T")
                row.names(IUPAC) = c("A", "G", "C", "T")
                return(as.character(IUPAC[l[1], l[2]]))
        }
        if (length(l) == 3) {
                if (!("A" %in% l)) {return("B")}
                if (!("C" %in% l)) {return("D")}
                if (!("G" %in% l)) {return("H")}
                if (!("T" %in% l)) {return("V")}
        }
}

for (i in files) {
        name_i = gsub("_lofreq.vcf", "", i)
        vcf_file = paste0(path, name_i, "_lofreq.vcf")
        cons_file = paste0(path, name_i, "_cons_gaps.fasta")
        cov_file = paste0(path, name_i, "_cov.list")
        
        vcf_data = read.table(vcf_file, quote="\"") %>%
                rename(CHROM = V1, POS = V2, ID = V3, REF = V4, ALT = V5, QUAL = V6, FILTER = V7, INFO = V8) %>%
                separate (INFO, INFO_lofreq, sep = ";", extra = "drop") %>%
                select(POS, REF, ALT, DP, AF) %>%
                mutate(DP = gsub("DP=" ,"", DP)) %>%
                mutate(AF = round(as.numeric(gsub("AF=" ,"", AF))*100,0)) %>%
                filter(AF >= threshold)
    
        seq_data = data.frame(cons = unlist(strsplit(readLines(cons_file)[2], ""))) %>%
                mutate(POS = seq.int(nrow(.)))
        
        cov_data = data.frame(COV = readLines(cov_file)) %>%
                mutate(POS = seq.int(nrow(.)))
        
        comb_data = full_join(seq_data, vcf_data, "POS") %>%
                full_join(cov_data, "POS") %>%
                mutate(REF = as.character(REF)) %>%
                mutate(cons = as.character(cons))
        
        if (!(all(comb_data$cons == comb_data$REF, na.rm = TRUE))) {next} ### Exit loop for this sample if alignment not correct
        
        comb_data = comb_data %>%
                rename(SEQ = cons) %>%
                select(POS, SEQ, ALT, AF, COV) %>%
                mutate(AF = ifelse(is.na(AF), 0, AF)) %>%
                mutate(ALT = as.character(ALT)) %>%
                mutate(ALT = ifelse(is.na(ALT), SEQ, ALT)) %>%
                mutate(WTS = apply(.[,c('SEQ', 'ALT')], 1, function(x) wobble(c(x['SEQ'], x['ALT'])))) #%>%
                #mutate(WTS = apply(.[,c('COV')], 1, function(x) ifelse(COV <= threshold, "-", WTS)))

        for (j in nrow(comb_data):2) {
                if (comb_data$POS[j] == comb_data$POS[j-1]) {
                        comb_data$WTS[j-1] = wobble(c(comb_data$SEQ[j], comb_data$SEQ[j-1], comb_data$ALT[j], comb_data$ALT[j-1]))
                        comb_data$AF[j-1] = paste(comb_data$AF[j], comb_data$AF[j-1], sep = "/")
                        comb_data$ALT[j-1] = paste(comb_data$ALT[j], comb_data$ALT[j-1], sep = "/")
                        comb_data = comb_data[-j, ]
                }
        }
        
        write.csv(comb_data, paste0(path, name_i,  "_", threshold, ".csv"))
        tr_w_seq = paste(comb_data$WTS, collapse = "")
        write.fasta(tr_w_seq, name_i, paste0(path, name_i, "_", threshold, "_wobble.fasta"), open = "w", nbchar = 60)
}
