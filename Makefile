ASSEMBLY_FILES = $(wildcard *-SOAPdenovo-Trans-assembly.fa)

all:
	@echo 'naw'

download:
	python scripts/download.py data/wanted_species.txt

clean: 
	rm -rf data/*-SOAPdeno-Trans-assembly.fa

%/: data/%-SOAPdenovo-Trans-assembly.fa
	makeblastdb -in $< -dbtype nucl -parse_seqids -out $@
	#mkdir -p $@-blastdb/
	#mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $@-blastdb/

.PHONY: all clean
.DELETE_ON_ERROR: