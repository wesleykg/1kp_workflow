all: $(patsubst data/%-SOAPdenovo-Trans-assembly.fa, data/%-blastdb, $(wildcard data/*-SOAPdenovo-Trans-assembly.fa))

data/*-SOAPdenovo-Trans-assembly.fa: data/wanted_species.txt
	cd scripts/ ; python download.py ../$?

data/%-blastdb: data/%-SOAPdenovo-Trans-assembly.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/PHIPA-genes/PHIPA_*.fasta: data/PHIPA-genes.fasta
	cd scripts/ ; python fasta_splitter.py $?
	mkdir data/PHIPA-genes/
	mv data/PHIPA_*.fasta data/PHIPA-genes/

cleantemp: 
	cd data/ ; rm -drf *-genes/

clean:
	cd data/ ; rm -drf *-blastdb/ *-genes/

cleanall: 
	cd data/ ; rm -drf *-SOAPdenovo-Trans-assembly.fa *-blastdb/ *-genes/

.PHONY: all clean cleantemp cleanall
.DELETE_ON_ERROR:
.PRECIOUS: data/%-SOAPdenovo-Trans-assembly.fa