all: $(patsubst data/%-SOAPdenovo-Trans-assembly.fa, data/%-blastdb, $(wildcard data/*-SOAPdenovo-Trans-assembly.fa))

download:
	cd scripts/;python download.py ../data/wanted_species.txt

clean:
	rm -rf data/*-blastdb/

cleanall: 
	rm -rf data/*-SOAPdeno-Trans-assembly.fa data/*-blastdb/

data/%-blastdb: data/%-SOAPdenovo-Trans-assembly.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

.PHONY: all clean download
.DELETE_ON_ERROR:
.PRECIOUS: data/%-SOAPdenovo-Trans-assembly.fa