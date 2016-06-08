blastdb: $(patsubst data/%-assembly.fa, data/%-blastdb, \
	$(wildcard data/*-assembly.fa))

search: $(patsubst data/$(genome)_%.fasta, data/%_*.xml, \
	$(wildcard data/*_*.fasta))

download:
	cd data/ ; python ../scripts/0_download.py wanted_species.txt

split:
	cd data/ ; python ../scripts/0_fasta-splitter.py $(genome)-genes.fasta

data/%-assembly-filtered.fa: data/%-assembly.fa data/%-stats.tsv
	cd data/ ; Rscript ../scripts/1_trs-filter.R $(notdir $^)

data/%-blastdb: data/%-assembly-filtered.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/%_*.xml: data/*_%.fasta
	cd data/ ; python ../scripts/2_search.py $(notdir $^)

cleantemp:
	cd data/ ; rm -drf *.xml blast-results.csv

clean:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml blast-results.csv

cleanall:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml blast-results.csv \
	*-assembly.fa *-stats.tsv

.PHONY: clean cleantemp cleanall download search split
.DELETE_ON_ERROR:
.PRECIOUS: data/%-assembly.fa data/%-stats.tsv data/%-assembly-filtered.fa