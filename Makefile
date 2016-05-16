all: $(patsubst data/%-SOAPdenovo-Trans-assembly.fa, data/%-blastdb,\
	 $(wildcard data/*-SOAPdenovo-Trans-assembly.fa))

download:
	cd data/ ; python ../scripts/download.py wanted_species.txt

data/%-blastdb: data/%-SOAPdenovo-Trans-assembly.fa
	mkdir -p $@
	cd data/ ; makeblastdb -in $(^F) -dbtype nucl -parse_seqids -out $*
	cd data/ ; mv $*.nhr $*.nin $*.nog $*.nsd $*.nsi $*.nsq $(@F)

data/PHIPA_*.fasta: data/PHIPA-genes.fasta
	cd data/ ; python ../scripts/fasta_splitter.py $(notdir $^)

data/PHIPA_%.txt: data/PHIPA_%.fasta data/*-blastdb/
	cd data/ ; blastn -query $^ -db $(word 2,$^) -evalue 0.001 -num_threads 4 \
	-out $(@F) -outfmt '6 qseqid sseqid pident qlen slen length evalue'

cleantemp: 
	cd data/ ; rm -drf *-genes/

clean:
	cd data/ ; rm -drf *-blastdb/ *-genes/

cleanall: 
	cd data/ ; rm -drf *-SOAPdenovo-Trans-assembly.fa *-blastdb/ *-genes/

.PHONY: all clean cleantemp cleanall download
.DELETE_ON_ERROR:
.PRECIOUS: data/%-SOAPdenovo-Trans-assembly.fa