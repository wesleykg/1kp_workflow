'''Usage: 2_search.py <query_gene>'''

# Modules
import os  # Manipulating filenames
from glob import glob  # Finding relevant databases to search
from Bio import SeqIO  # Reading the query sequence for blast
from Bio.Blast import NCBIXML  # Parsing blastn results
from Bio.Blast.Applications import NcbiblastnCommandline  # Running blastn
from Bio.Blast.Applications import NcbitblastxCommandline  # Running tblastx


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


def blastn(query, evalue, db, out):
    blastn_search = NcbiblastnCommandline(query=query, evalue=evalue,
                                          db=db, num_threads=2,
                                          out=blastn_xml_name, outfmt=5)
    blastn_search()  # initializes the search
    print 'Searching for', gene_name, 'in', db_name, 'using blastn'


def tblastx(query, evalue, db, out):
    tblastx_search = NcbitblastxCommandline(query=query, evalue=evalue,
                                            db=db, num_threads=2,
                                            out=tblastx_xml_name, outfmt=5)
    tblastx_search()  # initializes the search
    print 'Searching for', gene_name, 'in', db_name, 'using tblastx'

# There should be at least one blast database in
# the same folder with the name: 'ID-blastdb/ID'

# Reads in the query sequence for blastn and record the name and length of the
# gene
query_seq = SeqIO.read(query_file, 'fasta')
query_len = str(len(query_seq))
query_name = query_seq.id
gene_name = query_name.split('_')[1]

# Create results file if it doesn't exist and add headers to first row
# explaining the data
if os.path.exists(query_name + '_blast-results.csv') is False:
    with open(query_name + '_blast-results.csv', 'w') as results_file:
        results_file.write('query' + ',' + 'blast-db' + ',' + 'scaf' + ',' +
                           'query_len' + ',' + 'scaf_len' + ',' + 'align_len' +
                           ',' + 'e_val' + ',' + 'blast_ver' + ',' + 'seq' +
                           '\n')

# Loop through each database, and search for matches to the query sequence
# using blastn. Results for each search are written to an xml file.
for db in db_list:
    db_name = os.path.split(db)[1]
    db_id = db_name.split('-')[0]
    db_path = db + '/' + db_id

    blastn_xml_name = query_name + '_' + db_id + '_blastn' + '.xml'

    # Search using blastn, and save results to an xml file
    blastn(query_file, evalue=1e-20, db=db_path, out=blastn_xml_name)

    # blastn xml file read into memory here
    with open(blastn_xml_name, 'r') as blastn_xml:
        blastn_results = NCBIXML.read(blastn_xml)

    # Check if there were any hits found using blastn, and if not, do tblastx
    if not blastn_results.alignments:  # Syntax for checking if list is empty

        tblastx_xml_name = query_name + '_' + db_id + '_tblastx' + '.xml'

        tblastx(query_file, evalue=1e-20, db=db_path, out=tblastx_xml_name)

        # tblastx xml file read into memory here
        with open(tblastx_xml_name, 'r') as tblastx_xml:
                tblastx_results = NCBIXML.read(tblastx_xml)

                # Check if there were any hits found using tblastx
                if not tblastx_results.alignments:

                    with open(query_name + '_blast-results.csv', 'a') as \
                            out_results:
                        print 'No hits found for', gene_name, 'in', db_name
                        out_results.write(query_name + ',' + db_id + '\n')

        # Loop throuch each tblastx match for a single search and record the
        # name and length of the matching scaffold, alignment length, e-value,
        # and the DNA sequence of the matched hit
        for record in tblastx_results.alignments:
            scaf_len = str(record.length)
            scaf_name = record.accession
            for hsps in record.hsps:
                ali_len = str(hsps.align_length)
                e_val = str(hsps.expect)
                scaf_seq = hsps.sbjct

            # Append the results of a single tblastx match to
            # gene_blast-results.csv'
            with open(query_name + '_blast-results.csv', 'a') as out_results:
                out_results.write(query_name + ',' + db_id + ',' +
                                  scaf_name + ',' + query_len + ',' +
                                  scaf_len + ',' + ali_len + ',' + e_val +
                                  ',' + 'tblastx' + ',' + scaf_seq + '\n')

    # Loop throuch each match for a single blastn search and record the
    # name and length of the matching scaffold, alignment length, e-value,
    # and the DNA sequence of the matched hit
    for record in blastn_results.alignments:
        scaf_len = str(record.length)
        scaf_name = record.accession
        for hsps in record.hsps:
            ali_len = str(hsps.align_length)
            e_val = str(hsps.expect)
#            scaf_start_pos = str(hsps.sbjct_start)
#            scaf_end_pos = str(hsps.sbjct_end)
            scaf_seq = hsps.sbjct

            # Append the results of a single blastn match to
            # '_blast-results.csv'
            with open(query_name + '_blast-results.csv', 'a') as out_results:
                out_results.write(query_name + ',' + db_id + ',' +
                                  scaf_name + ',' + query_len + ',' +
                                  scaf_len + ',' + ali_len + ',' + e_val + ','
                                  'blastn' + ',' + scaf_seq + '\n')
