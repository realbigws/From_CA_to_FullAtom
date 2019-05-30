# Example for 'residue' objects

from modeller import *
from modeller.scripts import complete_pdb

def analyze_seq(description, seq):
    """Simple 'analysis' of a sequence of residues, from a model or alignment"""
    numcys = 0
    for res in seq:
        if res.pdb_name == 'CYS':
            numcys += 1
    print "%s contains %d residues, of which %d are CYS" \
          % (description, len(seq), numcys)

env = environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

mdl = complete_pdb(env, "1fas")

# 'mdl.residues' is a list of all residues in the model
print "1-letter code of 1st residue: ", mdl.residues[0].code
print "PDB name of residue '10': ", mdl.residues['10:'].pdb_name

# Get the first aligned sequence from a file
aln = alignment(env, file='../commands/toxin.ali')
firstseq = aln[0]

# Analyze all residues in the model, a subset, and all residues in the
# alignment sequence
analyze_seq("Model 1fas", mdl.residues)
analyze_seq("First 10 residues of 1fas", mdl.residue_range('1:', '10:'))
analyze_seq("Aligned sequence %s" % firstseq.code, firstseq.residues)
