# Example for: alignment.read_one()

from modeller import *
env = environ()

# Create an empty alignment
aln = alignment(env)

# Open the input alignment file, and get a handle to it:
input = modfile.File('toxin.ali', 'r')
# Same for the output file:
output = modfile.File('toxin-filter.ali', 'w')

# Read sequences one by one from the file handle in PIR format:
while aln.read_one(input, alignment_format='PIR'):
    print "Read code %s" % aln[0].code
    # Write only X-ray structures to the output file:
    if aln[0].prottyp == 'structureX':
        aln.write(output, alignment_format='FASTA')

# Explicitly close the files (not strictly necessary in this simple
# example, because they'll be closed at the end of the script anyway):
input.close()
output.close()
