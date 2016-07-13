blastdb: $(patsubst data/%-assembly.fa, data/%-blastdb, \
	$(wildcard data/*-assembly.fa))

search: $(patsubst data/$(genome)_%.fasta, data/$(genome)_%_*.xml, \
	$(wildcard data/*_*.fasta))

download:
	cd data/ ; python ../scripts/0_download.py wanted_species.txt

split:
	cd data/ ; python ../scripts/0_fasta_splitter.py $(genome).fasta

catenate:
	cd data/ ; cat *-assembly_cleaned.fasta > all_assembly_cleaned.fasta

data/%-assembly_cleaned.fasta: data/%-assembly.fa data/%-stats.tsv
	cd data/ ; Rscript ../scripts/1_trs_cleaner.R $(notdir $^)

data/%-blastdb: data/%-assembly_cleaned.fasta
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/$(genome)_%_*.xml: data/$(genome)_%.fasta
	cd data/ ; python ../scripts/2_search.py $(notdir $^)

cleantemp:
	cd data/ ; rm -drf *.xml *_blast-results.csv

clean:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml *_blast-results.csv

cleanall:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml *_blast-results.csv \
	*-assembly.fa *-assembly_cleaned.fa *-stats.tsv

.PHONY: clean cleantemp cleanall download search split
.DELETE_ON_ERROR:
.PRECIOUS: data/%-assembly.fa data/%-stats.tsv data/%-assembly_cleaned.fasta
