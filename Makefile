ASSEMBLY_FILES = $(wildcard *-SOAPdenovo-Trans-assembly.fa)

all:
	@echo 'naw'

clean: 
	rm -rf wanted_assembly_URLs.txt *.fa *.bz2

$(ASSEMBLY_FILES): data/wanted_species.txt
	python scripts/build_url.py data/wanted_species.txt

%/: $(ASSEMBLY_FILES)
	makeblastdb -in $< -dbtype nucl -parse_seqids -out $@
	#mkdir -p $@-blastdb/
	#mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $@-blastdb/

.PHONY: all clean
.DELETE_ON_ERROR: