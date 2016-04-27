'''Usage: retrieve.py <wanted_species>'''

from docopt import docopt # For command-line arguments
import requests

cmdln_args = docopt(__doc__) # Creates a dictionary of command-line arguments

wanted_species_file = cmdln_args.get('<wanted_species>')

prefix = "http://onekp.westgrid.ca/1kp-data/"
infix = "/assembly/"
suffix = "-SOAPdenovo-Trans-assembly.fa.bz2"

URL_list = []

with open(wanted_species_file, 'r') as species_file:
    for line in species_file:
        ID = line.split('-')[0]
        URL = prefix + line.rstrip() + infix + ID + suffix
        URL_list.append(URL)

for URL in URL_list:
    assembly = requests.get(URL, auth=('1kp-data', '1kp-rna1'))
    print assembly.text