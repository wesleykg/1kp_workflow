'''Usage: 2_search.py <query_gene>'''

# Modules
import os  # Manipulating filenames
from glob import glob  # Finding relevant databases to search
from Bio import SeqIO  # Reading the query sequence for blast
from Bio.Blast import NCBIXML  # Parsing blastn results
from Bio.Blast.Applications import NcbiblastnCommandline  # Running blastn


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
    query_file = cmdln_args.get('<query_gene>')
    db_list = glob(os.getcwd() + '/*-blastdb')
# Run interatively in an iPython console
if in_ipython() is True:
    query_file = '../data/PHYPA_accD.fasta'
    db_list = glob('../data/*-blastdb')

# There should be at least one blast database in
# the same folder with the name: 'ID-blastdb/ID'

# Reads in the query sequence for blastn and records the name and length of the
# gene
query_seq = SeqIO.read(query_file, 'fasta')
query_len = str(len(query_seq))
query_name = query_seq.id
gene_name = query_name.split('_')[1]

# Loop through each database, and search for matches to the query sequence
# using blastn. Results for each search are written as an xml file.
for db in db_list:
    db_name = os.path.split(db)[1]
    db_id = db_name.split('-')[0]
    db_path = db + '/' + db_id
    xml_name = query_name + '_' + db_id + '.xml'

    # Intialize the blastn module and set the parameters.
    blastn_search = NcbiblastnCommandline(query=query_file, evalue=0.001,
                                          db=db_path, num_threads=4,
                                          out=xml_name, outfmt=5)
    blastn_search()
    print 'Searching for', gene_name, 'in', db_name

    # Results from above blastn search are saved in an xml file and read into
    # memory here.
    with open(xml_name, 'r') as xml_results:
        results = NCBIXML.read(xml_results)

    # Loop throuch each match for a single blastn search and record the name
    # and length of the matching scaffold. The expect-value and the length of
    # the alignment are also recorded
    for record in results.alignments:
        scaf_len = str(record.length)
        scaf_name = record.accession
        for hsps in record.hsps:
            ali_len = str(hsps.align_length)
            e_val = str(hsps.expect)
            scaf_start_pos = str(hsps.sbjct_start)
            scaf_end_pos = str(hsps.sbjct_end)
            scaf_seq = hsps.sbjct

            # Append the results of a single match to 'blast-results.csv'
            with open('blast-results.csv', 'a') as out_results:
                out_results.write(query_name + ',' + query_len + ',' +
                                  scaf_name + ',' + scaf_len + ',' + ali_len +
                                  ',' + e_val + ',' + scaf_start_pos + ',' +
                                  scaf_end_pos + ',' + scaf_seq + '\n')
