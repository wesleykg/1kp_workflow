'''Usage: download.py <wanted_species>'''

from docopt import docopt #For command-line arguments
import requests #For downloading the assembly
from requests.auth import HTTPDigestAuth #For authenticating the connection
import bz2 #For decompressing the downloaded assembly

cmdln_args = docopt(__doc__) #Creates a dictionary of command-line arguments

#File should have a single 1kp id on each line, 
#following this form: ID-Genus_species
wanted_species_file = cmdln_args.get('<wanted_species>')

#Sets the standard parts of the URL to download from. The URL is produced as
#follows: base_url + ID-Genus_species + infix_url + ID + suffix_url
base_url = "http://onekp.westgrid.ca/1kp-data/"
infix_url = "/assembly/"
suffix_url = "-SOAPdenovo-Trans-assembly.fa.bz2"

outname_suffix = "-SOAPdenovo-Trans-assembly.fa"

#Open the list of wanted 1kp names and loop over each line, producing the ID 
#and species name and then pasting them together to create the URL. Requests
#downloads the assembly using Digest factor authentication. The assembly is
#then decompressed and saved using the original filename on the server.
with open(wanted_species_file, 'r') as species_file:
    for line in species_file:
        ID = line.split('-')[0]
        name_1kp = line.rstrip()
        URL = base_url + name_1kp + infix_url + ID + suffix_url
        print 'Downloading', line.rstrip()
        compressed_assembly = requests.get(URL, 
                                auth=HTTPDigestAuth('1kp-data', '1kp-rna1'))
        print 'Decompressing', line.rstrip()
        assembly = bz2.decompress(compressed_assembly.content)
        with open(ID + outname_suffix, 'w') as out_assembly:
            out_assembly.write(assembly)