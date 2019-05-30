# Example for: alignment.edit()

# Read an alignment, write it out in the 'PAP' format, with overhangs cut.

from modeller import *

log.level(1, 1, 1, 1, 0)
env = environ()
env.io.atom_files_directory = ['.', '../atom_files']

aln = alignment(env, file='overhang.ali', align_codes='all',
                alignment_format='PIR')

# Cut overhangs in the 1is4 sequence that are longer than 3 residues
# relative to the longest remaining entry in the alignment:
aln.edit(edit_align_codes='1is4', base_align_codes='rest',
         min_base_entries=1, overhang=3)
aln.write(file='overhang-1.pir', alignment_format='PIR')
aln.write(file='overhang-1.pap', alignment_format='PAP')
