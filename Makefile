blastdb: $(patsubst data/%-assembly.fa, data/%-blastdb, \
	$(wildcard data/*-assembly.fa))

search: $(patsubst data/$(genome)-%.fasta, data/$(genome)-%_blast-results.csv, \
	$(wildcard data/*-*.fasta))

align: $(patsubst data/$(genome)-%_blast-results.csv, \
	data/$(genome)-%_blast-alignment.fasta, \
	$(wildcard data/*-*_blast-results.csv))

download:
	cd data/ ; python ../scripts/0_download.py wanted_species.txt

split:
	cd data/ ; python ../scripts/0_fasta_splitter.py $(genome).fasta

catenate:
	cd data/ ; cat *-assembly_cleaned.fasta > all_assemblies_cleaned.fasta

data/%-assembly_cleaned.fasta: data/%-assembly.fa data/%-stats.tsv
	cd data/ ; Rscript ../scripts/1_trs_cleaner.R $(notdir $^)

data/%-blastdb: data/%-assembly_cleaned.fasta
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/$(genome)-%_blast-results.csv: data/$(genome)-%.fasta
	cd data/ ; python ../scripts/2_search.py $(notdir $^)

data/$(genome)-%_blast-alignment.fasta: data/$(genome)-%_blast-results.csv
	cd data/ ; python ../scripts/3_align.py $(notdir $^) all_assemblies_cleaned.fasta

cleantemp:
	cd data/ ; rm -drf *_blast-alignment.fasta

clean:
	cd data/ ; rm -drf *-*.fasta *-blastdb/ *_*.xml *_blast-results.csv \ 
	*_blast-alignment.fasta

cleanall:
	cd data/ ; rm -drf *-*.fasta *-blastdb/ *_*.xml *_blast-results.csv \
	*_blast-alignment.fasta *-assembly.fa *-assembly_cleaned.fa *-stats.tsv

.PHONY: align catenate clean cleantemp cleanall download search split
.DELETE_ON_ERROR:
.PRECIOUS: data/%-assembly.fa data/%-stats.tsv data/%-assembly_cleaned.fasta
