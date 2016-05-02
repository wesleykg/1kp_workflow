all: $(patsubst data/%-SOAPdenovo-Trans-assembly.fa, data/%-blastdb, $(wildcard data/*-SOAPdenovo-Trans-assembly.fa))

download:
	cd scripts/;python download.py ../data/wanted_species.txt

data/%-blastdb: data/%-SOAPdenovo-Trans-assembly.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/PHIPA-genes/PHIPA_%.fasta: data/PHIPA-genes.fasta
	cd scripts/;python fasta_splitter.py $?
	mkdir data/PHIPA-genes/
	mv data/PHIPA_*.fasta data/PHIPA-genes/

clean:
	rm -rf data/*-blastdb/ data/*_*.fasta

cleanall: 
	rm -rf data/*-SOAPdenovo-Trans-assembly.fa data/*-blastdb/ data/*_*.fasta

.PHONY: all clean cleanall download
.DELETE_ON_ERROR:
.PRECIOUS: data/%-SOAPdenovo-Trans-assembly.fa