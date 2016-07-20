'''Usage: 3_align.py <blast_results> <all_assemblies>'''

# Modules
import os  # Manipulating filenames
import pandas  # Reading in csv blast results
from glob import glob  # Finding assemblies to index
from Bio import SeqIO  # Indexing all scaffolds
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord  # Producing empty SeqRecords


# Check if running interactively in an iPython console, or in a script
# from the command-line
def in_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
# Run in a script from the command-line
if in_ipython() is False:
    from docopt import docopt  # Command-line argument handler
    cmdln_args = docopt(__doc__)
    blast_results_filename = cmdln_args.get('<blast_results>')
    assembly_list = glob(os.getcwd() + '/*-assembly_cleaned.fasta')
# Run interatively in an iPython console
if in_ipython() is True:
    blast_results_filename = '../data/PODLA-accD_blast-results.csv'
    assembly_list = glob('../data/*-assembly_cleaned.fasta')

query_name = os.path.split(blast_results_filename)[1]
query_name = query_name.split('_')[0]

# Index the concatenated file of all 1kp assemblies of interest
all_scaffolds = SeqIO.index_db('all_assemblies_index.idx',
                               filenames=assembly_list, format='fasta')

# Names of each column in the csv blast results file
table_header = 'query', 'blast_db', 'scaf', 'query_len', 'scaf_len', \
                'align_len', 'e_val', 'blast_ver', 'seq'

# Read in the csv blast results file
blast_results = pandas.read_csv(blast_results_filename, names=table_header)

# Retrieve results for blast searches that found a scaffold
hit_results = blast_results[blast_results.scaf != 'None found']

# Retrieve results for blast searches that didn't find a scaffold
missing_results = blast_results[blast_results.scaf == 'None found']

# Initalize list of wanted scaffolds to write to file
wanted_scaffold_seqs = []

# Loop through each scaffold name and search for the accompanying SeqRecord
# in the index of all scaffolds. Add each matched SeqRecord to the list of
# wanted scaffold sequences
wanted_scaffold_names = hit_results.scaf
for name in wanted_scaffold_names:
    name = name.rstrip()
    wanted_scaffold_seq = all_scaffolds[name]
    wanted_scaffold_seqs.append(wanted_scaffold_seq)

query_seq = SeqIO.read('../data/' + query_name + '.fasta', format='fasta')
wanted_scaffold_seqs.append(query_seq)

# Retrieve 1kp species IDs for blast searches with no scaffolds found. Loop
# through this list of ID's and create an empty SeqRecord for each ID
# missing a scaffold
missing_scaffold_names = missing_results.blast_db
for name in missing_scaffold_names:
    missing_seq_record = SeqRecord(Seq(''), id=name, description='')
    wanted_scaffold_seqs.append(missing_seq_record)

SeqIO.write(wanted_scaffold_seqs, query_name + '_blast-alignment.fasta',
            format='fasta')
