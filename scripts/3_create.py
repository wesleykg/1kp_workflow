'''Usage: 3_create.py <blast_results>'''

# Modules
import os  # Manipulating filenames
import pandas  # Reading in csv blast results
from glob import glob  # Finding assemblies to index
from Bio import SeqIO  # Indexing all scaffolds


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
    blast_results_filename = '../data/PODTO-psbM_blast-results.csv'
    assembly_list = glob('../data/*-assembly_cleaned.fasta')


# Check if an index of all scaffolds from all assemblies has been created
if os.path.exists('../data/all_assemblies_index.idx') is False:
    # Index all 1kp assemblies of interest
    all_scaffolds = SeqIO.index_db('../data/all_assemblies_index.idx',
                                   filenames=assembly_list, format='fasta')
else:
    # Read in the already existing index
    all_scaffolds = SeqIO.index_db('../data/all_assemblies_index.idx')

# Names of each column in the csv blast results file
table_header = 'query', 'blast_db', 'hit', 'query_len', 'hit_len', \
                'align_len', 'e_val', 'orientation', 'blast_ver'

# Read in the csv blast results file
blast_results = pandas.read_csv(blast_results_filename, names=table_header,
                                header=None, index_col=False)

# Subset results for blast searches that found a sense hit
sense_results = blast_results[blast_results.orientation == 'sense']

# Subset results for blast searches that found an antisense hit
antisense_results = blast_results[blast_results.orientation == 'antisense']

# Initalize list of wanted hits to write to file
wanted_hits = []

# Loop through each scaffold name and search for the accompanying SeqRecord
# in the index of all scaffolds. Add each matched SeqRecord to the list of
# wanted scaffold sequences
wanted_sense_names = sense_results.hit
for name in wanted_sense_names:
    name = name.rstrip()
    wanted_hit_seq = all_scaffolds[name]
    wanted_hits.append(wanted_hit_seq)

# Add the original query sequence to file
query_name = os.path.split(blast_results_filename)[1]
query_name = query_name.split('_')[0]
query_seq = SeqIO.read('../data/' + query_name + '.fasta', format='fasta')
wanted_hits.append(query_seq)

wanted_antisense_names = antisense_results.hit
for name in wanted_antisense_names:
    name = name.rstrip()
    wanted_hit_seq = all_scaffolds[name]
    wanted_hit_seq = wanted_hit_seq.reverse_complement(id=True, name=True,
                                                       description=True)
    wanted_hits.append(wanted_hit_seq)

SeqIO.write(wanted_hits, query_name + '_blast-unaligned.fasta',
            format='fasta')

# Retrieve 1kp species IDs for blast searches with no hits found. Convert the
# pandas dataframe into a list, then convert the list to a string
# Retrieve results for blast searches that didn't find a hit
missing_results = blast_results[blast_results.hit == 'None found']
missing_hit_names = missing_results.blast_db
missing_hit_names = missing_hit_names.values.tolist()
missing_hit_names = "\n".join(missing_hit_names)
with open(query_name + '_blast-missing.txt', 'w') as missing_hits:
    missing_hits.write(missing_hit_names)
