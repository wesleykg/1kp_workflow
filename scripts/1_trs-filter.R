"Usage: trs_cleaner.R <transcriptome> <transratestats>" -> doc

## Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

cmd_args <- docopt::docopt(doc) # Produce a vector of command line arguments

## Assign command line arguments to variables
transcriptome_file <- unlist(cmd_args[1]) 
transratestats_file <- unlist(cmd_args[2])

onekp_ID <- substr(transcriptome_file, start = 1, stop = 4)
clean_transcriptome_out <- paste0(onekp_ID, "-assembly-filtered", ".fa")

## Read in a single species Transcriptome (T) data
Tdat <- readDNAStringSet(filepath = transcriptome_file)

## Read in Transrate Stats (TRS) data
TRSdat <- read.delim(file = transratestats_file, stringsAsFactors = FALSE, 
                    skip = 44) # First 44 lines are a summary of the data
names(TRSdat)[1] <- 'names' # Simplify 1st column name

## Produce a list of names of good quality scaffolds
good_scaffolds_names <- TRSdat %>% filter(quality == "good") %>%
  select(names) %>% unlist

## Produce a vector that is an index of scaffold quality
scaffolds_index <- names(Tdat) %in% good_scaffolds_names

## Retain all good scaffolds; discard bad scaffolds
Cdat <- Tdat[scaffolds_index]

## Write good scaffolds to file:
writeXStringSet(Cdat, clean_transcriptome_out)
