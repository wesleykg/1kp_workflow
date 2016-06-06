blastdb: $(patsubst data/%-SOAPdenovo-Trans-assembly.fa, data/%-blastdb, \
	$(wildcard data/*-SOAPdenovo-Trans-assembly.fa))

search: $(patsubst data/$(genome)_%.fasta, data/%_*.xml, \
	$(wildcard data/*_*.fasta))

download:
	cd data/ ; python ../scripts/0_download.py wanted_species.txt

split:
	cd data/ ; python ../scripts/0_fasta-splitter.py $(genome)-genes.fasta

data/%-blastdb: data/%-SOAPdenovo-Trans-assembly.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/%_*.xml: data/*_%.fasta
	cd data/ ; python ../scripts/2_search.py $(notdir $^)

cleantemp:
	cd data/ ; rm -drf *_*.fasta

clean:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml blast-results.csv

cleanall:
	cd data/ ; rm -drf *_*.fasta *-blastdb/ *_*.xml blast-results.csv \
	*-SOAPdenovo-Trans-assembly.fa *-SOAPdenovo-Trans-Transrate-stats.tsv

.PHONY: clean cleantemp cleanall download search split
.DELETE_ON_ERROR:
.PRECIOUS: data/%-SOAPdenovo-Trans-assembly.fa \
	data/%-SOAPdenovo-Trans-Transrate-stats.tsv