# Example for: chain.filter(), chain.write()

# This will read a PDB file (segment), and write out all of its chains
# satisfying the listed conditions into separate alignment files in the
# PIR format.

from modeller import *

env = environ()
mdl = model(env, file='../atom_files/pdb1lzd.ent')
for c in mdl.chains:
    if c.filter(minimal_chain_length=30, minimal_resolution=2.0,
                minimal_stdres=30, chop_nonstd_termini=True,
                structure_types='structureN structureX'):
        filename = '1lzd%s.chn' % c.name
        print "Wrote out", filename
        atom_file, align_code = c.atom_file_and_code(filename)
        c.write(filename, atom_file, align_code, format='PIR',
                chop_nonstd_termini=True)
