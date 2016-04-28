'''Usage: retrieve.py <wanted_species>'''

#from docopt import docopt # For command-line arguments
import requests
from requests.auth import HTTPDigestAuth
import bz2

#cmdln_args = docopt(__doc__) # Creates a dictionary of command-line arguments

wanted_species_file = "../data/wanted_species.txt"# cmdln_args.get('<wanted_species>')

prefix = "http://onekp.westgrid.ca/1kp-data/"
infix = "/assembly/"
suffix = "-SOAPdenovo-Trans-assembly.fa.bz2"
outname_suffix = "-SOAPdeno-Trans-assembly.fa"
data_dir = '../data/'

with open(wanted_species_file, 'r') as species_file:
    for line in species_file:
        ID = line.split('-')[0]
        URL = prefix + line.rstrip() + infix + ID + suffix
        compressed_assembly = requests.get(URL, 
                                auth=HTTPDigestAuth('1kp-data', '1kp-rna1'))
        assembly = bz2.decompress(compressed_assembly.content)
        with open(data_dir + ID + outname_suffix, 'w') as out_assembly:
            out_assembly.write(assembly)