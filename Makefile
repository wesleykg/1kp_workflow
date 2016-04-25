all: $(patsubst %.fasta, %.txt, $(wildcard $(query_id)-plastid-genes/*.fasta))

clean: 
	rm -rf blast-results.tsv blast-results/*.txt

$(query_id)-plastid-genes/%.txt: $(query_id)-plastid-genes/%.fasta
	blastn -query $^ -db $(db_id)-transcriptome-db/$(db_id) -evalue 0.001 -num_threads 4 -out $@ -outfmt '6 qseqid sseqid pident qlen slen length evalue'
	cat $@ >> blast-results.tsv
	rm $@
	blastn -query $^ -db $(db_id)-transcriptome-db/$(db_id) -evalue 0.001 -num_threads 4 -out $@
	mv $@ blast-results/

.PHONY: all clean
.DELETE_ON_ERROR: