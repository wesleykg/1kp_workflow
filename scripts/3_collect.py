import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


all_assemblies_filename = '../data/all_assemblies_cleaned.fasta'
blast_results_filename = '../data/PHYPA_atpI_blast-results.csv'

all_scaffolds = SeqIO.index(all_assemblies_filename, format='fasta')

table_header = 'query', 'blast_db', 'scaf', 'query_len', 'scaf_len', \
                'align_len', 'e_val', 'blast_ver', 'seq'

wanted_scaffolds = []

blast_results = pandas.read_csv(blast_results_filename, names=table_header)

missing_results = blast_results[blast_results.scaf == 'None found']
missing_scaffold_names = missing_results.blast_db

for name in missing_scaffold_names:
    missing_seq = Seq('')
    missing_seq_record = SeqRecord(missing_seq, id=name, description='')
    wanted_scaffolds.append(missing_seq_record)

hit_results = blast_results[blast_results.scaf != 'None found']
wanted_scaffold_names = hit_results.scaf

for name in wanted_scaffold_names:
    name = name.rstrip()
    wanted_scaffold = all_scaffolds[name]
    wanted_scaffolds.append(wanted_scaffold)

SeqIO.write(wanted_scaffolds, 'wanted_scaffolds.fasta', format='fasta')
