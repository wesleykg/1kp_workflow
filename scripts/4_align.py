
from Bio.Align.Applications import MuscleCommandline  # Aligning sequences


# Align sequenes using muscle
def muscle_align(input, outname):
    muscle_cmdline = MuscleCommandline(input=input, out=outname)
    muscle_cmdline()


muscle_align(wanted_scaffold_seqs, outname=query_name +
             '_blast-alignment.fasta')
