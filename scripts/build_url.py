'''Usage: retrieve.py <wanted_species>'''

from docopt import docopt # For command-line arguments

cmdln_args = docopt(__doc__) # Creates a dictionary of command-line arguments

wanted_species_file = cmdln_args.get('<wanted_species>')

prefix = "http://onekp.westgrid.ca/1kp-data/"
middle = "/assembly/"
suffix = "-SOAPdenovo-Trans-assembly.fa.bz2"

out_URL_list = "wanted_assembly_URLs.txt"

URL_list = []

with open(wanted_species_file, 'r') as species_file:
    for line in species_file:
        ID = line.split('-')[0]
        URL = prefix + line.rstrip() + middle + ID + suffix
        URL_list.append(URL)

with open(out_URL_list, "w") as URLs_out:
	for URL in URL_list:
		URLs_out.write(URL + '\n')