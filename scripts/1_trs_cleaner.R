'Usage: 1_trs_cleaner.R <transcriptome> <transratestats>' -> doc

## Load packages
library(tools)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

cmd_args <- docopt::docopt(doc) # Produce a vector of command line arguments

## Assign command line arguments to variables
transcriptome_file <- unlist(cmd_args[1])
#transcriptome_file <- 'WOGB-assembly.fa'
transratestats_file <- unlist(cmd_args[2])
#transratestats_file <- 'WOGB-stats.tsv'

onekp_ID <- substr(transcriptome_file, start = 1, stop = 4)
filename <- file_path_sans_ext(transcriptome_file)
extension <- file_ext(transcriptome_file)
clean_transcriptome_out <- paste0(filename, '_cleaned', '.fasta')

## Read in a single species Transcriptome (T) data
if(extension == 'fa'){
	Tdat <- readDNAStringSet(filepath = transcriptome_file)
} else if (extension == 'faa'){
	Tdat <- readAAStringSet(filepath = transcriptome_file)
}

## Record the total number of scaffolds in the assembly
number_of_scaffolds <- length(Tdat)

## Read in Transrate Stats (TRS) data
TRSdat <- read.delim(file = transratestats_file, stringsAsFactors = FALSE, 
                    skip = 44) # First 44 lines are a summary of the data
names(TRSdat)[1] <- 'names' # Simplify 1st column name

## Produce a list of names of good quality scaffolds
good_scaffolds_names <- TRSdat %>% filter(quality == 'good') %>%
  select(names) %>% unlist

## Produce a vector that is an index of scaffold quality
scaffolds_index <- names(Tdat) %in% good_scaffolds_names

## Retain all good scaffolds; discard bad scaffolds
Cdat <- Tdat[scaffolds_index]

## Record the number of remaining scaffolds after cleaning out bad scaffolds
number_of_good_scaffolds <- length(Cdat)

## Record proportion of scaffolds removed to one decimal place
proportion_of_scaffolds_removed <- round(100 - 
                                    number_of_good_scaffolds/number_of_scaffolds 
                                    * 100, digits = 1)

print(paste0('Started with ', number_of_scaffolds, ' scaffolds in ', 
             transcriptome_file))
print(paste0('Ended with ', number_of_good_scaffolds, ' scaffolds in ', 
             clean_transcriptome_out))
print(paste0(proportion_of_scaffolds_removed, '% of scaffolds removed from ', 
             filename))

## Write good scaffolds to file:
writeXStringSet(Cdat, clean_transcriptome_out)
