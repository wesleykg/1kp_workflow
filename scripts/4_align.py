'''Usage: 4_align.py <unaligned_alignment>'''

# Modules
import os  # Manipulating filenames
from Bio.Align.Applications import MuscleCommandline  # Aligning sequences


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
    unaligned_alignment = cmdln_args.get('<unaligned_alignment>')
# Run interatively in an iPython console
if in_ipython() is True:
    unaligned_alignment = '../data/PHYPA-accD_blast-alignment.fasta'


# Retrieve filename for writing to file
alignment_name = os.path.split(unaligned_alignment)[1]
alignment_name = alignment_name.rpartition('_')[0]

muscle_align = MuscleCommandline(input=unaligned_alignment,
                                 out=alignment_name + '_aligned.fasta')
muscle_align()
