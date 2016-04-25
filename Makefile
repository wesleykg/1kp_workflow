all:
	echo 'naw'

clean: 
	rm -rf wanted_assembly_URLs.txt *.fa *.bz2

%-SOAPdenovo-Trans-assembly.fa: data/wanted_species.txt
	python scripts/build_url.py data/wanted_species.txt
	scripts/download.sh 

%/: %-SOAPdenovo-Trans-assembly.fa
	makeblastdb -in $< -dbtype nucl -parse_seqids -out $@
	#mkdir -p $@-blastdb/
	#mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $@-blastdb/

.PHONY: all clean
.DELETE_ON_ERROR: