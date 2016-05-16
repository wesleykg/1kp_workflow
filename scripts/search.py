'''Usage: search.py <query_gene> <db_location>'''

import os
from glob import glob
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
#from docopt import docopt

#cmdln_args = docopt(__doc__) #Creates a dictionary of command-line arguments

query_file = '../data/PHIPA_accD.fasta' #cmdln_args.get('<query_gene>')
db_list = glob('../data/*-blastdb') #cmdln_args.get('<query_gene>')

query_seq = SeqIO.read(query_file, 'fasta')
query_len = str(len(query_seq))
query_name = query_seq.id
gene_name = query_name.split('_')[1]

for db in db_list:
    db_name = os.path.split(db)[1]
    db_id = db_name.split('-')[0]
    db_path = db + '/' + db_id
    xml_name = gene_name + '_' + db_id + '.xml'
    blastn_search = NcbiblastnCommandline(query = query_file, evalue = 0.001,
                                      db = db_path, num_threads = 4,
                                      out = xml_name, outfmt = 5)
    blastn_search()
    print 'Searching for', gene_name, 'in', db_name

    with open(xml_name, 'r') as xml_results:
        results = NCBIXML.read(xml_results)
    
    for record in results.alignments:
        scaf_len = str(record.length)
        scaf_name = record.accession
        for hsps in record.hsps:
            ali_len =  str(hsps.align_length)
            e_val = str(hsps.expect)
            with open('blast-results.csv', 'a') as out_results:
                out_results.write(query_name + ',' + query_len + ',' +
                scaf_name + ',' + scaf_len + ',' + ali_len + ',' + e_val +
                '\n')