# Read in scaf_name column of *_blast-results.csv
# Produce index of ALL transcriptome file -> index fasta.index(filepath=glob, seqtype='DNA')
# Pull seqs from transcriptome file that match scaf_name from *_blast-results.csv

## Load packages
#library(tools)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

blast_results_header <- c('query', 'blast-db', 'scaf', 'query_len', 'scaf_len', 
													'align_len', 'e_val', 'blast_ver', 'seq')

blast_results <- read.csv(file = 'all_blast_results.csv', 
													stringsAsFactors = FALSE, 
													#na.strings = 'None found',
													col.names = blast_results_header)

wanted_scaffold_names <- blast_results %>% select(scaf) %>% 
  filter(scaf != 'None found') %>% unlist() %>% unname()

write(wanted_scaffold_names, file = 'wanted_scaffoldnames.txt')

names

all_scaffolds_index <- fasta.index(filepath = 'all_assemblies_cleaned.fasta')

wanted_scaffolds <- readDNAStringSet(all_scaffolds_index$desc[wanted_scaffold_names])

wanted_scaffolds <- all_scaffolds_index$desc[wanted_scaffold_names]

myvars <- c("v1", "v2", "v3")

#readDNAStringSet(all_assemblies_index)

