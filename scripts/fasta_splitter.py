'''Usage: fasta_splitter.py <fasta_file>'''

from Bio import SeqIO
import os

from docopt import docopt # For command-line arguments
cmdln_args = docopt(__doc__) # Creates a dictionary of command-line arguments
in_fasta_path = cmdln_args.get('<fasta_file>')

in_fasta = os.path.split(in_fasta_path)[1]

for record in list(SeqIO.parse(in_fasta, 'fasta')):
    filename = record.id + '.fasta'
    with open(filename, 'w') as out_fasta:
        SeqIO.write(record, out_fasta, 'fasta')