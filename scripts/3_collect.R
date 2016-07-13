# Read in scaf_name column of *_blast-results.csv
# Read in transcriptome file as index fasta.index(filepath=, seqtype='DNA')
# Pull seqs from transcriptome file that match scaf_name from *_blast-results.csv

read.csv(file = PHYPA_accD_blast-results.csv, stringsAsFactors = FALSE)

