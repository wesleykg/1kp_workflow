from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

query = 'PHIPA_accD.fasta'
db = 'RDOO-blastdb/RDOO'
xml_name = 'PHIPA_accD.xml'
data_dir = '../data/'

query_seq = SeqIO.read(data_dir + query, 'fasta')
query_len = str(len(query_seq))
query_name = query_seq.id

blastn_search = NcbiblastnCommandline(query = data_dir + query, 
                                      evalue = 0.001, db = data_dir + db, 
                                      num_threads = 4, 
                                      out = data_dir + xml_name, outfmt = 5)
blastn_search()

results = NCBIXML.read(open(data_dir + xml_name))

for record in results.alignments:
    scaf_len = str(record.length)
    scaf_name = record.title.rstrip('No definition line')
    for hsps in record.hsps:
        ali_len =  str(hsps.align_length)
        e_val = str(hsps.expect)

with open('blast-results.csv', 'a') as out_results:
    out_results.write(query_name + ',' + query_len + ',' + scaf_name + ',' +
    scaf_len + ',' + ali_len + ',' + e_val + '\n')