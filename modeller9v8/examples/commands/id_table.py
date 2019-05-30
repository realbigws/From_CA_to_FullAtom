# Example for: alignment.id_table(), alignment.compare_sequences(),
#              misc.principal_components(), misc.dendrogram()

# Pairwise sequence identity between sequences in the alignment.

from modeller import *

env = environ()
env.io.atom_files_directory = ['../atom_files']
# Read all entries in this alignment:
aln = alignment(env, file='toxin.ali')

# Access pairwise properties:
s1 = aln[0]
s2 = aln[1]
print "%s and %s have %d equivalences, and are %.2f%% identical" % \
      (s1, s2, s1.get_num_equiv(s2), s1.get_sequence_identity(s2))

# Calculate pairwise sequence identities:
aln.id_table(matrix_file='toxin_id.mat')

# Calculate pairwise sequence similarities:
mdl = model(env, file='2ctx', model_segment=('1:A', '71:A'))
aln.compare_sequences(mdl, rr_file='$(LIB)/as1.sim.mat', max_gaps_match=1,
                      matrix_file='toxin.mat', variability_file='toxin.var')
mdl.write(file='2ctx.var')

# Do principal components clustering using sequence similarities:
env.principal_components(matrix_file='toxin.mat', file='toxin.princ')

# Dendrogram in the log file:
env.dendrogram(matrix_file='toxin.mat', cluster_cut=-1.0)
