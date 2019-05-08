args <- commandArgs(TRUE)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(seqinr)

# set path to data (terminal slash is important!)
#path <- "/Volumes/data/Diagnostics/experiments/170815/"
path <- paste0(args[1],"/")

# minority variant threshold (%)
variant_threshold <- 15

# minimal coverage required (reads)
minimal_coverage <- 3

files <- list.files(path, pattern = "lofreq.vcf")

# grep the N_lofreq.vcf files with the highest iteration
files <- files[grep(paste0(max(str_sub(files, -12, -12)), "_lofreq.vcf"), files)]

call_wobbles <- function(nuc_list) {
    ### function to call ambiguous nucleotides (wobbles)
    ### input is a list of nucleotides (nuc_list, length may be greater than 2)

    nuc_list <- nuc_list %>%
        toupper() %>%
        unique() %>% .[!is.na(.)] # NAs in nuc_list are ignored (e.g. if position not covered, returns NA)

    for (i in 1:length(nuc_list)) {
        if (nchar(nuc_list[i]) > 1) {return("insertion")} # returns "insertion" if nuc_list not composed of single characters
        if (!(is.element(nuc_list[i], c("A", "G", "C", "T", "N")))) {return("unkown_nucleotide")} # returns "unkown_nucleotide" if characters other than A, G, C, T, N in nuc_list
    }

    if (length(nuc_list) == 0) {return(NA)}
    if ("N" %in% nuc_list | length(nuc_list) == 4) {return("N")} # returns "N" if at least one N or all four nucleotides included in nuc_list
    if (length(nuc_list) == 1) {return(nuc_list[1])} # one single nucleotide
    if (length(nuc_list) == 2) { # two different nucleotides
        IUPAC = data.frame(c("A", "R", "M", "W"), c("R", "G", "S", "K"), c("M", "S", "C", "Y"), c("W", "K", "Y", "T"))
        names(IUPAC) = c("A", "G", "C", "T")
        row.names(IUPAC) = c("A", "G", "C", "T")
        return(as.character(IUPAC[nuc_list[1], nuc_list[2]]))
    }
    if (length(nuc_list) == 3) { # three different nucleotides
        if (!("A" %in% nuc_list)) {return("B")}
        if (!("C" %in% nuc_list)) {return("D")}
        if (!("G" %in% nuc_list)) {return("H")}
        if (!("T" %in% nuc_list)) {return("V")}
    }
}

# Loop over all files
for (i in files) {

    name_i <- gsub("_lofreq.vcf", "", i)
    vcf_file <- paste0(path, name_i, "_lofreq.vcf")
    cons_file <- paste0(path, name_i, "_cons.fasta")
    depth_file <- paste0(path, name_i, ".depth")

    if (class(try(read.table(vcf_file))) == "try-error") {
        vcf_data = data.frame(POS = 1, REF = NA, ALT = NA, DP = NA, AF = NA) # if vcf file is empty
    } else {
      vcf_data = try(read_delim(vcf_file, "\t", escape_double = FALSE, col_names = FALSE,
                                col_types = cols(X4 = col_character(),
                                                 X5 = col_character()),
                                comment = "#",
                                trim_ws = TRUE)) %>%
        rename(CHROM = X1, POS = X2, ID = X3, REF = X4, ALT = X5, QUAL = X6, FILTER = X7, INFO = X8) %>%
        separate(INFO, c("DP", "AF", "SB", "DP4"), sep = ";", extra = "drop") %>%
        select(POS, REF, ALT, DP, AF) %>%
        mutate(DP = gsub("DP=" ,"", DP)) %>%
        mutate(AF = round(as.numeric(gsub("AF=" ,"", AF))*100,1)) %>%
        filter(AF >= variant_threshold) ### filter for AF >= variant_threshold
    }

    cons_data = data.frame(CONS = unlist(strsplit(readLines(cons_file)[-1], ""))) %>%
        mutate(POS = seq.int(nrow(.)))

    if (file.info(depth_file)$size == 0) {
        cov_data = data.frame(POS = 1, COV = NA) # if samtools.depth file is empty
    } else {
        cov_data = read_delim(depth_file, delim = "\t", col_names = FALSE, trim_ws = TRUE, col_types = "cii") %>%
            rename(POS = X2, COV = X3) %>%
            select(POS, COV)
    }

    comb_data = full_join(cons_data, vcf_data, "POS") %>%
        full_join(cov_data, "POS") %>%
        mutate(REF = toupper(as.character(REF))) %>%
        mutate(ALT = toupper(as.character(ALT))) %>%
        mutate(CONS = toupper(as.character(CONS)))

    ### filter for short amplicons
    #comb_data = comb_data %>% filter(POS >= 0 & POS <= 2000)

    ### call wobbles
    comb_data <- comb_data %>%
        mutate(REF = ifelse(is.na(REF), CONS,       ### use CONS when REF is NA
                            ifelse(REF == CONS, REF,       ### use REF (or CONS) if REF and CONS are identical
                                   ifelse(abs(AF - 50) <= 15, REF,  ### use REF although CONS is different when similar frequency
                                          "error")))) %>%                ### else "error"
        mutate(WTS = apply(.[,c('REF', 'ALT')], 1, function(x) call_wobbles(c(x['REF'], x['ALT'])))) %>%
        select(POS, REF, ALT, AF, COV, WTS)

    ### resolve duplicate lines in vcf_data for two variants
    for (j in nrow(comb_data):2) {
        if (comb_data$POS[j] == comb_data$POS[j - 1]) {
            comb_data$WTS[j - 1] = call_wobbles(c(comb_data$REF[j], comb_data$REF[j - 1], comb_data$ALT[j], comb_data$ALT[j - 1]))
            comb_data$AF[j - 1] = paste(comb_data$AF[j], comb_data$AF[j - 1], sep = "/")
            comb_data$ALT[j - 1] = paste(comb_data$ALT[j], comb_data$ALT[j - 1], sep = "/")
            comb_data = comb_data[-j, ]
        }
    }

    # apply minimal coverage threshold
    comb_data = comb_data %>%
        mutate(WTS = ifelse(is.na(COV), "-", ifelse(COV < minimal_coverage, "N", WTS)))

    # write output files
    write.csv(comb_data, paste0(path, name_i,  "_", variant_threshold, ".csv"))
    write.fasta(paste(comb_data$WTS, collapse = ""), name_i, paste0(path, name_i, "_", variant_threshold, "_WTS.fasta"), open = "w", nbchar = 60)
}
