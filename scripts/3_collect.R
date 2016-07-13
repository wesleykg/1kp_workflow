# Read in scaf_name column of *_blast-results.csv
# Produce index of ALL transcriptome file -> index fasta.index(filepath=glob, seqtype='DNA')
# Pull seqs from transcriptome file that match scaf_name from *_blast-results.csv

## Load packages
#library(tools)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

blast_results_header <- c('query', 'blast-db', 'scaf', 'query_len', 'scaf_len', 
													'align_len', 'e_val', 'blast_ver', 'seq')

blast_results <- read.csv(file = 'PHYPA_all_blast-results.csv', 
													stringsAsFactors = FALSE, 
													col.names = blast_results_header)

wanted_scaffold_names <- blast_results %>% select(scaf)

all_scaffolds_index <- fasta.index(filepath = 'all-assembly_cleaned.fasta')

wanted_scaffolds <- all_scaffolds_index[wanted_scaffold_names]

#readDNAStringSet(all_assemblies_index)

