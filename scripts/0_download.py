'''Usage: 0_download.py <wanted_species>'''

# Modules
import requests  # Downloading the assemblies and stats files
from requests.auth import HTTPDigestAuth  # Authenticate the connection
import bz2  # Decompressing the downloaded assembly
from gzip import GzipFile  # Decompressing the downloaded stats file
from StringIO import StringIO  # Reading the compressed stats file


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
    from docopt import docopt  # Command-line arguments handler
    cmdln_args = docopt(__doc__)
    wanted_species_file = cmdln_args.get('<wanted_species>')
# Run interatively in an iPython console
if in_ipython() is True:
    wanted_species_file = '../data/wanted_species.txt'


# Sets the standard parts of the assembly_url to download from. The
# assembly_url is produced as
# follows: base_assembly_url + ID-Genus_species + infix_assembly_url + ID +
# assembly_suffix
base_assembly_url = 'http://206.12.96.204/1kp-data/'
infix_assembly_url = '/assembly/'
assembly_suffix = '-SOAPdenovo-Trans-assembly.fa.bz2'
stats_suffix = '-SOAPdenovo-Trans-Transrate-stats.tsv.gz'

outname_assembly = '-assembly.fa'
outname_stats = '-stats.tsv'

# Open the list of wanted 1kp names and loop over each line, producing the ID
# and species name and then pasting them together to create the necessary URLs.
# Both files are downloaded using Digest-factor authentication. The files are
# then decompressed and saved using a truncated version of the server filename.
with open(wanted_species_file, 'r') as species_file:
    for line in species_file:
        ID = line.split('-')[0]
        name_1kp = line.rstrip()
        species_url = base_assembly_url + name_1kp + infix_assembly_url + ID
        assembly_url = species_url + assembly_suffix
        stats_url = species_url + stats_suffix
        print 'Downloading', line.rstrip(), 'assembly'
        compressed_assembly = requests.get(assembly_url, auth=HTTPDigestAuth
                                           ('1kp-data', '1kp-rna1'))
        print 'Decompressing', line.rstrip(), 'assembly'
        assembly = bz2.decompress(compressed_assembly.content)
        with open(ID + outname_assembly, 'w') as out_assembly:
            out_assembly.write(assembly)
        print 'Downloading', line.rstrip(), 'statistics file'
        compressed_stats = requests.get(stats_url, auth=HTTPDigestAuth
                                        ('1kp-data', '1kp-rna1'))
        print 'Decompressing', line.rstrip(), 'statistics file'
        stats = GzipFile(fileobj=StringIO(compressed_stats.content))
        with open(ID + outname_stats, 'w') as out_stats:
            out_stats.write(stats.read())
